import scipy.sparse as sp
import numpy as np
from collections import defaultdict
import scipy.integrate as integrate

class String:
    def __init__(self, config: dict) -> None:
        self.config = config
        self.__get_M_sh()
        self.__get_D_h()
        self.__get_M_qh()


    def __get_M_sh(self) -> None:
        Ns = self.config["Ns"]
        ls = self.config["ls"]
        rho_s = self.config["rho_s"]
        diag_values = np.full(Ns - 1, ls * rho_s / (Ns - 1))
        self.M_sh = sp.diags(diag_values, offsets=0, format='coo')

    def __get_D_h(self) -> None:
        Ns = self.config["Ns"]
        T = self.config["T"]
        row = np.arange(Ns - 1)
        col = np.arange(Ns - 1)
        vals = np.full(Ns - 1, -T)
        row_extra = np.arange(Ns - 1)
        col_extra = np.arange(1, Ns)
        vals_extra = np.full(Ns - 1, T)
        rows = np.concatenate((row, row_extra))
        cols = np.concatenate((col, col_extra))
        values = np.concatenate((vals, vals_extra))
        self.D_h = sp.coo_matrix((values, (rows, cols)), shape=(Ns - 1, Ns))

    def __h_t(self, t) -> float:
        assert t >= 0.0, "t must be larger or equal to zero"
        t1 = self.config["t1"]
        t2 = self.config["t2"]
        if 0.0 <= t < t1:
            return 1.0 - np.cos(np.pi * t / t1)
        elif t < t2:
            return 1.0 + np.cos(np.pi * (t - t1) / (t2 - t1))
        else:
            return 0.0

    def get_f_sh(self, t) -> None:
        x0 = self.config["x0"]
        delta_s = self.config["delta_s"]
        ls = self.config["ls"]
        Ns = self.config["Ns"]
        integrand = lambda x: np.exp(-((x - x0) / delta_s) ** 2)
        int_0_ls, error = integrate.quad(integrand, 0.0, ls)
        inv_int_0_ls = 1.0 / int_0_ls
        ht = self.__h_t(t)
        ret = np.array([integrate.quad(integrand, i * ls / (Ns - 1), (i + 1) * ls / (Ns - 1))[0] for i in range(Ns - 1)])
        ret *= ht * inv_int_0_ls
        self.f_sh = ret

    def __get_M_qh(self) -> None:
        ls = self.config["ls"]
        Ns = self.config["Ns"]
        T = self.config["T"]
        Q0 = 1 / (3 * T) * (ls / (Ns - 1))
        Q1 = 1 / (6 * T) * (ls / (Ns - 1))
        rows = []
        cols = []
        values = []
        for i in range(Ns - 1):
            rows.extend([i, i + 1, i, i + 1])
            cols.extend([i, i + 1, i + 1, i])
            values.extend([Q0, Q0, Q1, Q1])
        self.M_qh = sp.coo_matrix((values, (rows, cols)), shape=(Ns, Ns))

from grid import Grid3d
import scipy.sparse as sp
from collections import defaultdict
import scipy.integrate as integrate
from utils import create_sprase_mat

class Air(Grid3d):

    def __init__(self, config):
        self.config = config
        grid_min = np.array(config["grid min"])
        grid_delta = np.array(config["grid resolution"]) * config["grid l"]
        grid_max = grid_min + grid_delta
        super().__init__(grid_min, grid_max, config["grid resolution"])
        self.__get_M_pah()
        self.__get_Gh_transpose()
        self.__get_M_ah()
    

    def __get_M_pah(self) -> None:
        diag_values = np.full(self.get_num_grid() , self.config["grid l"] ** 3 / (self.config["rho_a"] * self.config["c_a"] ** 2))
        self.M_pah = sp.diags(diag_values, offsets=0, format="coo")

    def __get_Gh_transpose(self)->None:
        mat_data:defaultdict[tuple[int, int], float] = defaultdict(lambda : 0.0)

        h = self.config["grid l"]
        mat_int = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)
        for j in range(self.get_num_grid()):
            global_idxs = self.get_idx1d_base_func_idxs(j)
            for i in range(6):
                mat_data[(j, global_idxs[i])] += h**2 * mat_int[i]
        
        self.Gh_transpose = create_sprase_mat(mat_data, (self.get_num_grid(), self.get_num_base_func()))

    def __get_M_ah(self) -> None:
        mat_int:np.ndarray[np.ndarray[float]] = np.array([
            [2.0, 1.0, 0.0, 0.0, 0.0, 0.0], 
            [1.0, 2.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 2.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 2.0, 1.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 2.0]
        ]) * self.config["grid l"] ** 3 / 6 * self.config["rho_a"]

        mat_data:defaultdict[tuple[int, int], float] = defaultdict(lambda : 0.0)
        
        for grid_idx in range(self.get_num_grid()):
            global_idxs = self.get_idx1d_base_func_idxs(grid_idx)
            for j in range(6):
                for i in range(6):
                    mat_data[(global_idxs[i], global_idxs[j])] += mat_int[i, j]

        self.M_ah = create_sprase_mat(mat_data, (self.get_num_base_func(), self.get_num_base_func()))


from utils import map_vertices, map_edge, make_edge, get_boundary_edges, reorder_tris, screen_tris
from vis_utils import show_vector
import open3d as o3d

class SoundBoard:
    def __init__(
        self,
        vertices: np.ndarray[tuple[float, float, float]],
        tris: np.ndarray[tuple[int, int, int]],
        config: dict
    ) -> None:
        print("SoundBoard Info:")
        self.vertices = vertices
        inner_tris: np.ndarray[tuple[int, int, int]] = screen_tris(
            self.vertices, tris, config["plate vertical axis"], config["plate thd"])
        self.vertices_map = map_vertices(inner_tris)
        self.num_vertices = len(self.vertices_map)
        self.edge_map = map_edge(inner_tris, self.num_vertices)
        self.num_edges = len(self.edge_map)
        self.num_centers = len(inner_tris)
        self.config = config
        self.data = np.array([(Idx_tri_0, Idx_tri_1, Idx_tri_2, self.vertices_map[Idx_tri_0], self.vertices_map[Idx_tri_1], self.vertices_map[Idx_tri_2],
                               self.edge_map[make_edge(Idx_tri_0, Idx_tri_1)], self.edge_map[make_edge(
                                   Idx_tri_1, Idx_tri_2)], self.edge_map[make_edge(Idx_tri_0, Idx_tri_2)],
                               self.num_vertices + self.num_edges + i) for i, (Idx_tri_0, Idx_tri_1, Idx_tri_2) in enumerate(inner_tris)])
        self.gamma_f = reorder_tris(inner_tris, get_boundary_edges(vertices, inner_tris, config["plate x axis"], config["plate x min"], config[
                                    "plate x max"], config["plate y axis"], config["plate y min"], config["plate y max"], True), config["gamma_f_ref_dir"])
        self.gamma_0 = reorder_tris(inner_tris, get_boundary_edges(vertices, inner_tris, config["plate x axis"], config["plate x min"], config[
                                    "plate x max"], config["plate y axis"], config["plate y min"], config["plate y max"], False), config["gamma_0_ref_dir"])

        print(f"# of vertices: {self.num_vertices}")
        print(f"# of edges: {self.num_edges}")
        print(f"# of centriods: {self.num_centers}")
        print(f"# of total points: {self.total_points()}")
        print(f"# of edges on gamma_f: {len(self.gamma_f)}")
        print(f"# of edges on gamma_0: {len(self.gamma_0)}")
        self.__get_J_h()
        self.__get_H_h()
        self.__get_M_Mh()
        self.__get_M_ph()
        self.__get_B_w_h_transpose()

    def total_points(self) -> int:
        return self.num_vertices + self.num_edges + self.num_centers

    def get_gamma_f_edges(self, color: tuple[float, float, float] = (0.2, 0.2, 0.7)) -> list[o3d.geometry]:
        ret = []
        for idx_start, idx_end, tri_idx in self.gamma_f:
            ret.append(show_vector(
                self.vertices[idx_start], self.vertices[idx_end], color))
        return ret

    def get_gamma_0_edges(self, color: tuple[float, float, float] = (0.2, 0.2, 0.7)) -> list[o3d.geometry]:
        ret = []
        for idx_start, idx_end, tri_idx in self.gamma_0:
            ret.append(show_vector(
                self.vertices[idx_start], self.vertices[idx_end], color))
        return ret

    def __get_M_ph(self) -> None:
        M_int = [
            [13 / 240, 29 / 1440, 29 / 1440, 1 / 20, 1 / 45, 1 / 20, 3 / 56],
            [29 / 1440, 13 / 240, 29 / 1440, 1 / 20, 1 / 20, 1 / 45, 3 / 56],
            [29 / 1440, 29 / 1440, 13 / 240, 1 / 45, 1 / 20, 1 / 20, 3 / 56],
            [1 / 20, 1 / 20, 1 / 45, 4 / 45, 2 / 45, 2 / 45, 2 / 45, 3 / 35],
            [1 / 45, 1 / 20, 1 / 20, 2 / 45, 4 / 45, 2 / 45, 3 / 35],
            [1 / 20, 1 / 45, 1 / 20, 2 / 45, 2 / 45, 4 / 45, 3 / 35],
            [3 / 56, 3 / 56, 3 / 56, 3 / 35, 3 / 35, 3 / 35, 81 / 560]
        ]

        self.M_ph = sp.dok_matrix((self.total_points(), self.total_points()))
        a = self.config["a"]
        rho_p = self.config["rho_p"]

        for loc_i0, loc_i1, loc_i2, *global_idxs in self.data:
            p0 = self.vertices[loc_i0]
            p1 = self.vertices[loc_i1]
            p2 = self.vertices[loc_i2]

            vec_u = p1 - p0
            vec_v = p2 - p0
            vec_u[1] = 0.0
            vec_v[1] = 0.0
            u_cross_v = np.cross(vec_u, vec_v)
            coef = a * rho_p * np.linalg.norm(u_cross_v)

            for i in range(7):
                for j in range(7):
                    self.M_ph[global_idxs[j], global_idxs[i]] += coef * M_int[i][j]

    def __get_J_h(self) -> None:
        x0 = self.config["x0"]
        y0 = self.config["y0"]
        Ns = self.config["Ns"]
        T = self.config["T"]
        
        self.J_h = sp.dok_matrix((self.total_points(), Ns))

        f0 = lambda u, v: (1.0 - (u + v) / 2) * (1. - u - v)
        f1 = lambda u, v: (1.0 + u) * u / 2.0
        f2 = lambda u, v: (1.0 + v) * v / 2.0
        f3 = lambda u, v: 4.0 * u * (1.0 - u - v)
        f4 = lambda u, v: 4 * u * v
        f5 = lambda u, v: 4 * v * (1.0 - u - v)
        f6 = lambda u, v: 27 * u * v * (1.0 - u - v)

        f = [f0, f1, f2, f3, f4, f5, f6]

        cnt = 0
        for tri_idx, (loc_i0, loc_i1, loc_i2, *global_idxs) in enumerate(self.data):
            p0 = self.vertices[loc_i0]
            p1 = self.vertices[loc_i1]
            p2 = self.vertices[loc_i2]

            vec_target = np.array([x0, 0, y0]) - p0
            vec_u = p1 - p0
            vec_v = p2 - p0
            vec_target[1] = 0.0
            vec_u[1] = 0.0
            vec_v[1] = 0.0
            
            u = vec_v[2] * vec_target[0] - vec_v[0] * vec_target[2]
            v = vec_u[0] * vec_target[2] - vec_u[2] * vec_target[0]
            det = vec_u[0] * vec_v[2] - vec_u[2] * vec_v[0]
            u /= det
            v /= det

            if u >= 0.0 and v >= 0.0 and u + v <= 1.0:
                cnt += 1
                for i in range(7):
                    self.J_h[(global_idxs[i], Ns-1)] += f[i](u, v) * T

        if cnt != 1:
            print(cnt)

        assert cnt == 1, "only one satisfied triangle"

    def __get_M_Mh(self) -> None:
        M_int = [
            [13 / 240, 29 / 1440, 29 / 1440, 1 / 20, 1 / 45, 1 / 20, 3 / 56],
            [29 / 1440, 13 / 240, 29 / 1440, 1 / 20, 1 / 20, 1 / 45, 3 / 56],
            [29 / 1440, 29 / 1440, 13 / 240, 1 / 45, 1 / 20, 1 / 20, 3 / 56],
            [1 / 20, 1 / 20, 1 / 45, 4 / 45, 2 / 45, 2 / 45, 3 / 35],
            [1 / 45, 1 / 20, 1 / 20, 2 / 45, 4 / 45, 2 / 45, 3 / 35],
            [1 / 20, 1 / 45, 1 / 20, 2 / 45, 2 / 45, 4 / 45, 3 / 35],
            [3 / 56, 3 / 56, 3 / 56, 3 / 35, 3 / 35, 3 / 35, 81 / 560]
        ]

        D1 = self.config["D1"]
        D2 = self.config["D2"]
        D3 = self.config["D3"]
        D4 = self.config["D4"]

        inv_det = 1 / (D1 * D3 - D2 * D2 / 4)

        inv_C = [
            [D3 * inv_det, -D2 * inv_det / 2,   0],
            [-D2 * inv_det / 2, D1 * inv_det,   0],
            [0,                0, 4 / D4]
        ]

        data: dict[tuple[int, int], float] = defaultdict(lambda: 0.0)
        self.M_Mh = sp.dok_matrix((self.total_points() * 3, self.total_points() * 3))

        for loc_i0, loc_i1, loc_i2, *global_idx in self.data:
            loc_p0 = self.vertices[loc_i0]
            loc_p1 = self.vertices[loc_i1]
            loc_p2 = self.vertices[loc_i2]

            vec_u = loc_p1 - loc_p0
            vec_v = loc_p2 - loc_p0
            vec_u[1] = 0.0
            vec_v[1] = 0.0
            norm_cross_uv = np.linalg.norm(np.cross(vec_u, vec_v))

            for i in range(7):
                for j in range(i, 7):
                    int_v = M_int[i][j] * norm_cross_uv
                    for ii in range(3):
                        for jj in range(3):
                            b_ii_jj = inv_C[ii][jj]
                            self.M_Mh[(3 * global_idx[i] + ii, 3 * global_idx[j] + jj)] += int_v * b_ii_jj


    def __get_H_h(self) -> None:
        M_int_uu = [
            [3 / 8, -1 / 3, 0, -1 / 6, -1 / 2, 1 / 2, -9 / 40],
            [-1 / 3, 3 / 8, 0, -1 / 6, 1 / 2, -1 / 2, -9 / 40],
            [0, 0, 0, 0, 0, 0, 0],
            [-1 / 6, -1 / 6, 0, 4 / 3, 0, 0, 9 / 5],
            [-1 / 2, 1 / 2, 0, 0, 4 / 3, -4 / 3, 0],
            [1 / 2, -1 / 2, 0, 0, -4 / 3, 4 / 3, 0],
            [-9 / 40, -9 / 40, 0, 9 / 5, 0, 0, 81 / 20]
        ]

        M_int_uv = [
            [3 / 8, 0, -1 / 3, 1 / 2, -1 / 2, -1 / 6, -9 / 40],
            [-1 / 3, 0, 1 / 3, -2 / 3, 2 / 3, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [-1 / 6, 0, 0, 2 / 3, -2 / 3, 2 / 3, 9 / 10],
            [-1 / 2, 0, 2 / 3, -2 / 3, 2 / 3, -2 / 3, -9 / 10],
            [1 / 2, 0, -2 / 3, 2 / 3, -2 / 3, 2 / 3, 9 / 10],
            [-9 / 40, 0, 0, 9 / 10, -9 / 10, 9 / 10, 81 / 40]
        ]

        M_int_vv = [
            [3 / 8, 0, -1 / 3, 1 / 2, -1 / 2, -1 / 6, -9 / 40],
            [0, 0, 0, 0, 0, 0, 0],
            [-1 / 3, 0, 3 / 8, -1 / 2, 1 / 2, -1 / 6, -9 / 40],
            [1 / 2, 0, -1 / 2, 4 / 3, -4 / 3, 0, 0],
            [-1 / 2, 0, 1 / 2, -4 / 3, 4 / 3, 0, 0],
            [-1 / 6, 0, -1 / 6, 0, 0, 4 / 3, 9 / 5],
            [-9 / 40, 0, -9 / 40, 0, 0, 9 / 5, 81 / 20]
        ]

        base_func_map: dict[tuple[int, int], tuple[tuple[int, int, int], tuple[int, int, int]]] = {
            (0, 1): ([0, 1, 2], [0, 1, 2]),
            (1, 0): ([1, 0, 2], [0, 2, 1]),
            (0, 2): ([0, 2, 1], [2, 1, 0]),
            (2, 0): ([2, 0, 1], [2, 0, 1]),
            (1, 2): ([1, 2, 0], [1, 2, 0]),
            (2, 1): ([2, 1, 0], [1, 0, 2])
        }

        M_int_u = [
            [-1 / 2, 1 / 3, 0, 2 / 3, 0, 0, 0],
            [-1 / 3, 1 / 2, 0, -2 / 3, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [-2 / 3, 2 / 3, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0]
        ]

        self.H_h = sp.dok_matrix((self.total_points() * 3, self.total_points()))

        # plate

        for loc_i0, loc_i1, loc_i2, *global_idx in self.data:

            loc_p0 = self.vertices[loc_i0]
            loc_p1 = self.vertices[loc_i1]
            loc_p2 = self.vertices[loc_i2]

            vec_u = loc_p1 - loc_p0
            vec_v = loc_p2 - loc_p0
            vec_u[1] = 0.0
            vec_v[1] = 0.0
            norm_corss_uv = np.linalg.norm(np.cross(vec_u, vec_v))
            Q = np.array([
                [vec_u[0], vec_v[0]],
                [vec_u[2], vec_v[2]]
            ])
            inv_Q = np.linalg.inv(Q)

            for i in range(7):
                for j in range(7):
                    self.H_h[(3 * global_idx[j] + 0, global_idx[i])] += norm_corss_uv * (M_int_uu[i][j] * inv_Q[0, 0] * inv_Q[0, 0] + (
                        M_int_uv[i][j] + M_int_uv[j][i]) * inv_Q[0, 0] * inv_Q[1, 0] + M_int_vv[i][j] * inv_Q[1, 0] * inv_Q[1, 0])
                    self.H_h[(3 * global_idx[j] + 1, global_idx[i])] += norm_corss_uv * (M_int_uu[i][j] * inv_Q[0, 1] * inv_Q[0, 1] + (
                        M_int_uv[i][j] + M_int_uv[j][i]) * inv_Q[0, 1] * inv_Q[1, 1] + M_int_vv[i][j] * inv_Q[1, 1] * inv_Q[1, 1])
                    self.H_h[(3 * global_idx[j] + 2, global_idx[i])] += norm_corss_uv * (2.0 * M_int_uu[i][j] * inv_Q[0, 0] * inv_Q[0, 1] + (M_int_uv[i]
                                                                                                                                         [j] + M_int_uv[j][i]) * (inv_Q[0, 0] * inv_Q[1, 1] + inv_Q[0, 1] * inv_Q[1, 0]) + 2.0 * M_int_vv[i][j] * inv_Q[1, 0] * inv_Q[1, 1])

        # gamma_f
        for loc_i0, loc_i1, tri_idx in self.gamma_f:
            tri_data = self.data[tri_idx]
            loc_idxs, global_idx_012, global_idx_345 = list(
                tri_data[:3]), np.array(tri_data[3:6]), np.array(tri_data[6:9])
            transform_012, transform_345 = base_func_map[(
                loc_idxs.index(loc_i0), loc_idxs.index(loc_i1))]
            global_idx_012 = global_idx_012[transform_012]
            global_idx_345 = global_idx_345[transform_345]

            vec_tao = self.vertices[loc_i1] - self.vertices[loc_i0]
            vec_tao[1] = 0.0
            inv_len_vec_tao = 1. / np.linalg.norm(vec_tao)

            vec_n = np.cross(vec_tao, [0.0, 1.0, 0.0])
            vec_n /= np.linalg.norm(vec_n)

            global_idx = list(tri_data[3:])
            tri_data[3:6] = list(global_idx_012)
            tri_data[6:9] = list(global_idx_345)

            for j in range(7):
                for i in [0, 1, 3]:
                    self.H_h[(3 * global_idx[j] + 0, global_idx[i])] += inv_len_vec_tao * vec_n[0] * vec_tao[0] * M_int_u[i][j]
                    self.H_h[(3 * global_idx[j] + 1, global_idx[i])] += inv_len_vec_tao * vec_n[2] * vec_tao[2] * M_int_u[i][j]
                    self.H_h[(3 * global_idx[j] + 2, global_idx[i])] += inv_len_vec_tao * (vec_n[0] * vec_tao[2] + vec_n[2] * vec_tao[0]) * M_int_u[i][j]

    def __get_B_w_h_transpose(self):
        M_int_gf: list[list[float]] = [
            [1 / 15, 7 / 240, 7 / 240, 1 / 15, 1 / 30, 1 / 15, 3 / 40],
            [23 / 240, 23 / 240, 7 / 120, 2 / 15, 1 / 10, 1 / 10, 3 / 20],
            [23 / 240, 7 / 120, 23 / 240, 1 / 10, 1 / 10, 2 / 15, 3 / 20]
        ]

        self.B_w_h_transpose = sp.dok_matrix((self.total_points(), len(self.vertices)))

        for loc_i0, loc_i1, loc_i2, *global_idxs in self.data:
            p0 = self.vertices[loc_i0]
            p1 = self.vertices[loc_i1]
            p2 = self.vertices[loc_i2]

            vec_u = p1 - p0
            vec_v = p2 - p0
            vec_u[1] = 0.0
            vec_v[1] = 0.0

            u_cross_v = np.cross(vec_u, vec_v)
            len_u_cross_v = np.linalg.norm(u_cross_v)

            for i in range(7):
                for j in range(3):
                    self.B_w_h_transpose[global_idxs[i], global_idxs[j]] += len_u_cross_v * M_int_gf[j][i]


if __name__ == '__main__':
    import open3d as o3d
    import numpy as np
    mesh = o3d.io.read_triangle_mesh("sound_board.stl")

    from constants import constants
    vertices = np.asanyarray(mesh.vertices)
    triangles = np.asanyarray(mesh.triangles)

    vertices, inverse_indices = np.unique(vertices, axis=0, return_inverse=True)
    triangles = inverse_indices[triangles]

    from models import SoundBoard
    soundBoard = SoundBoard(vertices, triangles, constants)
