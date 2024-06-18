import matplotlib.pyplot as plt
import numpy as np
import math
import open3d as o3d




class Grid3d : 
    def __init__(self, xmin:float, xmax:float, ymin:float, ymax:float, zmin:float, zmax:float, x_res:int, y_res:int, z_res:int) -> None:
        self.min = np.array([xmin, ymin, zmin])
        self.max = np.array([xmax, ymax, zmax])
        assert(x_res >= 2 and y_res >= 2 and z_res >= 2)
        self.res = np.array([x_res, y_res, z_res])
        
    def get_delta(self, dim:int)->float:
        assert(0 <= dim < 3)
        return (self.max[dim] - self.min[dim]) / (self.res[dim] - 1)

    def get_deltas(self)->list[float, float, float]:

        return (self.max - self.min) / (self.res - 1)

    def get_1d_idx(self, pos:float, dim:int)->tuple[int, int]:
        assert(0 <= dim < 3)
        assert(self.min[dim] <= pos <= self.max[dim]) 

        delta = self.get_delta(dim)
        ceil = math.ceil(pos / delta)
        floor = math.floor(pos / delta)
        return (ceil, floor)
    
    def get_3d_idx(self, pos3d)->tuple[tuple[int, int], tuple[int, int] , tuple[int, int]]:
        return (self.get_1d_idx(pos3d[0], 0), self.get_1d_idx(pos3d[1], 1), self.get_1d_idx(pos3d[2], 2))
    
    def get_pos1d(self, idx:int, dim:int)->float:
        assert(0 <= dim < 3)
        assert(0 <= idx < self.res[dim])
        return self.min[dim] + idx * self.get_delta(dim)
        
    
    def get_pos3d(self, idx3d)->tuple[float, float, float]:
        return (self.get_pos1d(idx3d[0], 0), self.get_pos1d(idx3d[1], 1), self.get_pos1d(idx3d[2], 2))

    def vis_all_points(self)->None:

        points = np.zeros([self.res[0] * self.res[1] * self.res[2], 3])

        delta_x, delta_y, delta_z = self.get_delta(0), self.get_delta(1), self.get_delta(2)

        for z_idx in range(self.res[2]):
            z_pos = self.min[2] + delta_z * z_idx
            for y_idx in range(self.res[1]):
                y_pos = self.min[1] + delta_y * y_idx
                for x_idx in range(self.res[0]):
                    x_pos = self.min[0] + delta_x * x_idx
                    global_idx = x_idx + self.res[0] * (y_idx + self.res[1] * z_idx)
                    points[global_idx, 0] = x_pos
                    points[global_idx, 1] = y_pos
                    points[global_idx, 2] = z_pos

        point_cloud = o3d.geometry.PointCloud()
        point_cloud.points = o3d.utility.Vector3dVector(points)
        o3d.visualization.draw_geometries([point_cloud])

if __name__ == "__main__":
    grid3d = Grid3d(-0.4, 0.4, -0.2, 0.2, -0.05, 0.05, 17, 9, 3)
    grid3d.vis_all_points()
    pass

