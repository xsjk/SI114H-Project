import matplotlib.pyplot as plt
import numpy as np
import math
import open3d as o3d

class Grid3d:
    def __init__(self, 
                min:tuple[float, float, float], 
                max:tuple[float, float, float], 
                resolution:tuple[int, int, int]
                ) -> None:

        '''@resolution:在长、宽、高上grid的数量'''
        self.min:np.ndarray = np.array(min)
        self.max:np.ndarray = np.array(max)
        self.resolution = np.array(resolution)


    def __get_delta(self)->np.ndarray:
        delta = self.max - self.min
        delta /= self.resolution
        return delta

    def __idx_check(self, idx3d:tuple[int, int, int])->None:
        assert(0 <= idx3d[0] < self.resolution[0])
        assert(0 <= idx3d[1] < self.resolution[1])
        assert(0 <= idx3d[2] < self.resolution[2])

    def get_point_idx3d(self, pos:tuple[float, float, float])->tuple[int, int, int]:
        pos = np.array(pos)
        idx3d = (pos - self.min) / self.__get_delta()
        idx3d = np.floor(idx3d)
        idx3d = np.int32(idx3d)
        self.__idx_check(idx3d)
        return idx3d

    def get_point_idx1d(self, pos:tuple[float, float, float])->tuple[int]:
        idx3d = self.get_point_idx3d(pos)
        return idx3d[0] + self.resolution[0] * (idx3d[1] + self.resolution[1] * idx3d[2])

    def idx1d_to_idx3d(self, idx1d:int)->tuple[int, int, int]:
        x_idx = int(idx1d % self.resolution[0])
        y_idx = int((idx1d / self.resolution[0]) % self.resolution[1])
        z_idx = int((idx1d / self.resolution[0]) / self.resolution[2])
        return (x_idx, y_idx, z_idx)


    
    def get_idx3d_base_func_idx(self, idx3d:tuple[int, int, int])->tuple[int, int, int, int, int, int]:
        '''@return 
            ret[0]:xmin base function idx
            ret[1]:xmax base function idx
            ret[2]:ymin base function idx
            ret[3]:ymax base function idx
            ret[4]:zmin base function idx
            ret[5]:zmax base function idx
        '''
        x_idx, y_idx, z_idx = idx3d

        x_min_base = x_idx * self.resolution[1] * self.resolution[2]
        x_max_base = (x_idx + 1) * self.resolution[1] * self.resolution[2]
        x_offset = y_idx * self.resolution[2] + z_idx

        y_min_base = y_idx * self.resolution[0] * self.resolution[2]
        y_max_base = (y_idx + 1) * self.resolution[0] * self.resolution[2]
        y_offset = z_idx * self.resolution[0] + x_idx

        z_min_base = z_idx * self.resolution[0] * self.resolution[1]
        z_max_base = (z_idx + 1) * self.resolution[0] * self.resolution[1]
        z_offset = x_idx * self.resolution[1] + y_idx

        #local idx to global
        y_offset += self.resolution[1] * self.resolution[2] * (self.resolution[0] + 1)
        z_offset += self.resolution[1] * self.resolution[2] * (self.resolution[0] + 1)
        z_offset += self.resolution[0] * self.resolution[2] * (self.resolution[1] + 1)

        return (x_min_base + x_offset, x_max_base + x_offset, y_min_base + y_offset, y_max_base + y_offset, z_min_base + z_offset, z_max_base + z_offset)


    def get_idx1d_base_func_idxs(self, idx1d:int)->tuple[int, int, int, int, int, int, int]:
        return self.get_idx3d_base_func_idx(self.idx1d_to_idx3d(idx1d))

    def get_point_base_func_idxs(self, pos:tuple[float, float, float])->tuple[int, int, int, int, int, int]:
        return self.get_idx3d_base_func_idx(self.get_point_idx3d(pos))
    
    def get_num_base_func(self)->int:
        return 3 * self.resolution[0] * self.resolution[1] * self.resolution[2] + self.resolution[0] * self.resolution[1] + self.resolution[0] * self.resolution[2] + self.resolution[1] * self.resolution[2]
    
    def get_num_grid(self)->int:
        return self.resolution[0] * self.resolution[1] * self.resolution[2]
    
    def get_min_coord(self, pos:tuple[float, float, float])->tuple[float, float, float]:
        pos = np.array(pos)
        return self.min + np.floor((pos - self.min) / self.__get_delta()) * self.__get_delta()


if __name__ == "__main__":
    grid3d = Grid3d((0.0,0.0,0.0), (1.0, 1.0, 1.0), (1, 1, 1))
    print(grid3d.get_idx1d_base_func_idxs(0))



    

