from grid import Grid3d
from enum import Enum
import numpy as np
import open3d as o3d

class Air:

    class PointType(Enum):
        FREE = 1
        SOUND_BOARD_NEIGHBOR = 2
        RIGID_NEIGHBOR = 4

    PointColor:dict = {
        PointType.FREE : np.array([0.1, 0.1, 0.1], dtype = np.float32), 
        PointType.SOUND_BOARD_NEIGHBOR : np.array([0.9, 0.1, 0.1], dtype = np.float32), 
        PointType.RIGID_NEIGHBOR : np.array([0.1, 0.9, 0.1], dtype = np.float32)  
    }

    @staticmethod
    def is_free_point(pointType:PointType)->bool:
        return bool(pointType.value & Air.PointType.FREE.value)
    
    @staticmethod
    def is_sound_board_neighbor(pointType:PointType)->bool:
        return bool(pointType.value & Air.PointType.SOUND_BOARD_NEIGHBOR.value)

    @staticmethod
    def is_rigid_neighbor(pointType:PointType)->bool:
        return bool(pointType.value & Air.PointType.RIGID_NEIGHBOR.value)
    
    def __init__(self, xmin:float, xmax:float, ymin:float, ymax:float, zmin:float, zmax:float, x_res:int, y_res:int, z_res:int) -> None:
        self.grid = Grid3d(xmin, xmax, ymin, ymax, zmin, zmax, x_res , y_res, z_res)
        self.types:dict[Air.PointType, list[tuple[int, int, int]]] = dict()
        for type in Air.PointType:
            self.types[type] = []

    def set_point_type(self, point3d, type:PointType)->None:
        (xs, ys, zs) = self.grid.get_3d_idx(point3d)
        for z in zs:
            for y in ys:
                for x in xs:
                    self.types[type].append((x, y, z))
    
    def set_triangle_points(self, p0:list[float, float, float], p1:list[float, float, float], p2:list[float, float, float], type:PointType)->None:
        v0:list[float, float, float] = p1 - p0
        v1:list[float, float, float] = p2 - p0
        deltas = self.grid.get_deltas()
        delta_u:float = np.min(deltas / (np.abs(v0) + 1e-10))
        delta_v:float = np.min(deltas / (np.abs(v1) + 1e-10))
        v = 0.0
        while(v <= 1.0):
            vec_v = v * v1
            u = 0.0
            while(u <= 1.0 - v):
                vec_u = u * v0
                point = p0 + vec_v + vec_u
                self.set_point_type(point, type)
                u += delta_u
            v += delta_v

    def set_triangles_points(self, vertices:list[list[float, float, float]], indices:list[list[int, int, int]], type)->None:
        for (i0, i1, i2) in indices:
            self.set_triangle_points(vertices[i0], vertices[i1], vertices[i2], type)

    def set_points_type(self, points3d, type:PointType):
        for point3d in points3d:
            self.set_point_type(point3d, type)
    
    def get_point_clouds(self, type:PointType)->o3d.geometry.PointCloud:
        assert(type in Air.PointType)
        point_cloud = o3d.geometry.PointCloud()
        point_cloud.points = o3d.utility.Vector3dVector(np.array(self.types[type]))
        point_cloud.colors = o3d.utility.Vector3dVector(np.array([Air.PointColor[type]] * len(self.types[type])))
        return point_cloud
    
    




