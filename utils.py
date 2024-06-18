import numpy as np
from collections import defaultdict
import scipy.sparse as sp

def make_edge(i0: int, i1: int) -> tuple[int, int]:
    return (i0, i1) if i0 < i1 else (i1, i0)


def screen_tris(
    vertices: np.ndarray[tuple[float, float, float]],
    tris: np.ndarray[tuple[int, int, int]],
    dim: int,
    thd: float
) -> np.ndarray[tuple[int, int, int]]:

    assert (0 <= dim <= 2)
    index = np.array(range(len(tris)))
    return tris[(vertices[tris[index, 0], dim] < thd) & (vertices[tris[index, 1], dim] < thd) & (vertices[tris[index, 2], dim] < thd)]


def map_vertices(tris: np.ndarray[tuple[int, int, int]], offset: int = 0) -> dict[int, int]:
    ret = dict()

    for i0, i1, i2 in tris:
        if i0 not in ret:
            ret[i0] = offset + len(ret)
        if i1 not in ret:
            ret[i1] = offset + len(ret)
        if i2 not in ret:
            ret[i2] = offset + len(ret)

    return ret


def map_edge(
    tris: np.ndarray[tuple[int, int, int]],
    offset: int
) -> dict[tuple[int, int], int]:

    ret = dict()

    for i0, i1, i2 in tris:
        e0, e1, e2 = make_edge(i0, i1), make_edge(i1, i2), make_edge(i2, i0)
        if e0 not in ret:
            ret[e0] = offset + len(ret)
        if e1 not in ret:
            ret[e1] = offset + len(ret)
        if e2 not in ret:
            ret[e2] = offset + len(ret)

    return ret


def reverse_map(d: dict[int, int]) -> dict[int, int]:
    ret = dict()
    for i in d:
        ret[d[i]] = i

    return ret


def dict_create_or_add(d: dict, key: tuple[int, int], v: float) -> None:
    if key not in d:
        d[key] = v
    else:
        d[key] += v


def get_coef(v0, v1, p) -> tuple[float, float]:
    u = np.dot(v0, p) / np.linalg.norm(v0)
    v = np.dot(v1, p) / np.linalg.norm(v1)

    if u >= 0 and v >= 0 and (u + v) <= 1:
        return (u, v)
    else:
        return False


def get_boundary_edges(
        vertices: np.ndarray[tuple[float, float, float]],
        tris: np.ndarray[tuple[float, float, float]],
        axis0: int,
        min0: float,
        max0: float,
        axis1: int,
        min1: float,
        max1: float,
        is_inner: bool) -> list[tuple[int, int, int]]:
    # @return <edge_vertice_idx0, edge_vertice_idx1, tri_idx>

    edge_map: dict[tuple[int, int], list[int]] = defaultdict(lambda: [])
    # edge_map.key <edge_vertice_idx0, edge_vertice_idx1> such that edge_vertice_idx0 < edge_vertice_idx1
    # edge_map.value <list[tri_idx]>
    def judge_func(cord3d): return (min0 <= cord3d[axis0] <= max0 and min1 <= cord3d[axis1] <= max1) if is_inner else (
        cord3d[axis0] < min0 or cord3d[axis0] > max0 or cord3d[axis1] < min1 or cord3d[axis1] > max1)

    for tri_idx, (i0, i1, i2) in enumerate(tris):
        edge_map[make_edge(i0, i1)].append(tri_idx)
        edge_map[make_edge(i1, i2)].append(tri_idx)
        edge_map[make_edge(i0, i2)].append(tri_idx)

    ret = []
    for idx0, idx1 in edge_map:
        if len(tmp := edge_map[(idx0, idx1)]) == 1 and judge_func(vertices[idx0]) and judge_func(vertices[idx1]):
            ret.append((idx0, idx1, tmp[0]))
    return ret


def adapt_one_tri(
    tris: np.ndarray[tuple[int, int, int]],
    idx: int,
    target_i0: int,
    target_i1: int
) -> None:
    assert target_i0 in tris[idx], f"{target_i0} must in tris[{idx}]"
    assert target_i1 in tris[idx], f"{target_i1} must in tris[{idx}]"
    i2 = None
    for i in tris[idx]:
        if i != target_i0 and i != target_i1:
            i2 = i
            break
    tris[idx][0] = target_i0
    tris[idx][1] = target_i1
    tris[idx][2] = i2


def reorder_tris(
        tris: np.ndarray[tuple[int, int, int]],
        boundary: list[tuple[int, int, int]],
        ref_dir: bool = True) -> list[tuple[int, int, int]]:
    # @return list[<edge_start, edge_end, tri_idx>]
    edge_map: dict[int, list[tuple[int, int]]] = defaultdict(lambda: [])
    # edge_map.key <edge_vertice_idx>
    # edge_map.value <list[<edge_vertice_neighbor, tri_idx>]>

    for i0, i1, tri_idx in boundary:
        edge_map[i0].append((i1, tri_idx))
        edge_map[i1].append((i0, tri_idx))

    vis: set = set()
    ret: list[int] = []
    head = None
    next = None
    if ref_dir:
        head = boundary[0][0]
        next = boundary[0][1]
    else:
        head = boundary[0][1]
        next = boundary[0][0]
    vis.add(head)
    ret.append((head, next, boundary[0][2]))

    while next not in vis:
        (i0, tri0), (i1, tri1) = edge_map[next]
        vis.add(next)
        if i0 not in vis:
            ret.append((next, i0, tri0))
            next = i0
        if i1 not in vis:
            ret.append((next, i1, tri1))
            next = i1

    (i0, tri0), (i1, tri1) = edge_map[next]
    if i0 == head:
        ret.append((next, head, tri0))
    if i1 == head:
        ret.append((next, head, tri1))

    return ret



def compute_normals(vertices, triangles, reference_point):
    """Compute the (external) normals of the mesh.

    Args:
        vertices (array-like, 3d points): The vertices of the mesh.
        triangles (array-like, 3-tulple of indexes): The triangles of the mesh.
        reference_point (array-like, one 3d point): The point to which the normals will be directed away from.
    
    Return:
        normals (array-like, 3d vectors): The normals of the mesh.
    """
    # Please use more numpy functions to speed up the computation.
    vertices_num = len(vertices)
    triangles_num = len(triangles)
    normals = np.zeros((triangles_num, 3))
    for i in range(triangles_num):
        p1, p2, p3 = vertices[triangles[i]].reshape(3, 3)
        normal = np.cross(p2 - p1, p3 - p1)
        assert np.linalg.norm(normal) > 0
        normal /= np.linalg.norm(normal)
        if np.dot(normal, reference_point - p1) < 0:
            normal *= -1
        normals[i] = normal
    return normals

def check_normals(vertices, triangles, normals, reference_point) -> bool:
    """Check whether the normals are correct (external and directed away from the reference point).

    Args:
        vertices (array-like): The vertices of the mesh.
        triangles (array-like, 3-tuples): The triangles of the mesh.
        normals (array-like, 3d vectors): The normals of the mesh.
        reference_point (3d point): The point to which the normals should be directed away from.
    Return:
        bool: Whether the normals are correct.
    """
    for i in range(len(triangles)):
        p1, p2, p3 = vertices[triangles[i]].reshape(3, 3)
        normal = normals[i]
        edge1 = p2 - p1
        edge2 = p3 - p1
        if np.abs(np.dot(normal, edge1)) > 1e-6 or np.abs(np.dot(normal, edge2)) > 1e-6:
            return False
        if np.dot(normal, reference_point - p1) < 0:
            return False
    return True