import open3d as o3d
import numpy as np
from enum import Enum
from time import time

class Guitar:

    class PointType(Enum):
        sigma = 0
        omega = 1
        gamma_f = 2
        gamma_0 = 4
    
    @staticmethod
    def in_omega(type:PointType)->bool:
        return bool(type.value & Guitar.PointType.omega.value)
    
    @staticmethod
    def in_gamma_f(type:PointType)->bool:
        return bool(type.value & Guitar.PointType.gamma_f.value)
    
    @staticmethod
    def in_gamma_0(type:PointType)->bool:
        return bool(type.value & Guitar.PointType.gamma_0.value)
    
    @staticmethod
    def in_sigma(type:PointType)->bool:
        return bool(type.value == Guitar.PointType.sigma.value)
            


    def __init__(self, path:str) -> None:
        print("Loading Guitar Begin...")
        t_start = time()
        mesh = o3d.io.read_triangle_mesh(path)
        self.vertices = np.asanyarray(mesh.vertices)
        self.triangles = np.asanyarray(mesh.triangles)
        self.vertices_type = [Guitar.PointType.sigma for _ in range(len(self.vertices))]
        self.__remove_repeat()
        t_end = time()
        print("Loading Guitar Finish")
        print(f"# of vertices {len(self.vertices)}.")
        print(f"# of triangles {len(self.triangles)}.")
        print(f"loading time {t_end - t_start} sec.")

    @staticmethod
    def sort_index(i0:int, i1:int, i2:int)->tuple[int, int, int]:
        if i0 < i1 < i2:
            return (i0, i1, i2)
        elif i0 < i2 < i1:
            return (i0, i2, i1)
        elif i1 < i0 < i2:
            return (i1, i0, i2)
        elif i1 < i2 < i0:
            return (i1, i2, i0)
        elif i2 < i0 < i1:
            return (i2, i0, i1)
        else:
            return (i2, i1, i0)

    def __remove_repeat(self)->None:
        point_to_idx:dict[tuple[float, float, float], int] = dict()

        for (x, y, z) in self.vertices:
            if (x, y, z) not in point_to_idx:
                point_to_idx[(x, y, z)] = len(point_to_idx)
        
    
        new_vertices = np.array([list(key) for key in point_to_idx.keys()])
        new_triangles = np.array(
            list(set([
                Guitar.sort_index(
                    point_to_idx[tuple(self.vertices[i0])], 
                    point_to_idx[tuple(self.vertices[i1])], 
                    point_to_idx[tuple(self.vertices[i2])]
                )
                for (i0, i1, i2) in self.triangles
            ]))
        )
        self.vertices = new_vertices
        self.triangles = new_triangles

    def __get_omega(self, thd:float = -0.0464)->None:
        p0 = self.triangles.T[0]
        p1 = self.triangles.T[1]
        p2 = self.triangles.T[2]

        boolean_idx = (self.vertices[p0, 1] < thd) & (self.vertices[p1, 1] < thd) & (self.vertices[p2, 1] < thd)



        pass


if __name__ == "__main__":
    guitar = Guitar("sound_board.stl")
    pass