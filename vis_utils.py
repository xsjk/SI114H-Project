import open3d as o3d
import numpy as np
X_AXIS = o3d.geometry.TriangleMesh.create_arrow(
    cylinder_radius=0.0005, cone_radius=0.001, cylinder_height=0.008, cone_height=0.002
)
X_AXIS.compute_vertex_normals()
X_AXIS.paint_uniform_color([1.0, 0.0, 0.0])
X_AXIS.rotate(o3d.geometry.get_rotation_matrix_from_axis_angle(
    np.array([0.0, np.pi / 2, 0.0]).reshape(3, 1)), center=(0, 0, 0))
Y_AXIS = o3d.geometry.TriangleMesh.create_arrow(
    cylinder_radius=0.0005, cone_radius=0.001, cylinder_height=0.008, cone_height=0.002
)
Y_AXIS.compute_vertex_normals()
Y_AXIS.paint_uniform_color([0.0, 1.0, 0.0])
Y_AXIS.rotate(o3d.geometry.get_rotation_matrix_from_axis_angle(
    np.array([np.pi * 1.5, 0.0, 0.0]).reshape(3, 1)), center=(0, 0, 0))
Z_AXIS = o3d.geometry.TriangleMesh.create_arrow(
    cylinder_radius=0.0005, cone_radius=0.001, cylinder_height=0.008, cone_height=0.002
)
Z_AXIS.compute_vertex_normals()
Z_AXIS.paint_uniform_color([0.0, 0.0, 1.0])


def show_vector(start: tuple[float, float, float], end: tuple[float, float, float], color: tuple[float, float, float]) -> o3d.geometry:
    start = np.array(start)
    end = np.array(end)
    dir = end - start
    l_dir = np.linalg.norm(dir)
    e_dir = dir / l_dir

    cylinder_height = l_dir * 0.6
    cone_height = l_dir * 0.4
    cylinder_radius = cylinder_height / 16
    cone_radius = cone_height / 2
    ret = o3d.geometry.TriangleMesh.create_arrow(
        cylinder_radius=cylinder_radius, cone_radius=cone_radius, cylinder_height=cylinder_height, cone_height=cone_height
    )
    ret.compute_vertex_normals()
    ret.paint_uniform_color(color)
    rot_axis = np.cross(np.array([0.0, 0.0, 1.0]), e_dir)
    rot_axis /= np.linalg.norm(rot_axis)
    cos_theta = np.dot(np.array([0.0, 0.0, 1.0]), e_dir)  # since unit vector
    theta = np.arccos(cos_theta)
    theta += 2 * np.pi
    ret.rotate(o3d.geometry.get_rotation_matrix_from_axis_angle(
        rot_axis.reshape(3, 1) * theta), center=(0, 0, 0))
    ret.translate(start)
    return ret