
import numpy as np
import cupy as cp
import scipy.sparse as sp
import cupyx.scipy.sparse.linalg

from models import String, SoundBoard, Air, Gamma
from sparse_utils import cpu_to_gpu

class PhysicalSimulation:
    def __init__(self, constants, string: String, soundboard: SoundBoard, air: Air, gamma: Gamma):

        cpu_to_gpu(constants)

        self.string = string
        self.soundboard = soundboard
        self.air = air
        self.gamma = gamma

        self.constants = constants
        self.time_step = constants['delta_t']
        
        self.variables = {
            'v_s': cp.zeros((constants['Ns']-1, 1)),  # Piecewise constant on the string
            'q': cp.zeros((constants['Ns'], 1)),       # Piecewise linear on the string

            # Variables on the soundboard
            'v_p': cp.zeros((soundboard.num_vertices + soundboard.num_edges + soundboard.num_centers, 1)),
            'a_p': cp.zeros((soundboard.num_vertices + soundboard.num_edges + soundboard.num_centers, 1)), # the time derivative of v_p
            'M_p': cp.zeros(((soundboard.num_vertices + soundboard.num_edges + soundboard.num_centers) * 3, 1)),

            # Variables in the air
            'v_a': cp.zeros((3 * (air.resolution[0] * air.resolution[1] * air.resolution[2]) + 
                                  air.resolution[0] * air.resolution[1] + 
                                  air.resolution[1] * air.resolution[2] + 
                                  air.resolution[0] * air.resolution[2], 1)),
            'p_a': cp.zeros((air.resolution[0] * air.resolution[1] * air.resolution[2], 1)),

            # Variables on the guitar body
            'lambda': cp.zeros((len(gamma.vertices), 1))
        }

        # TODO set initial states

        self.history = {k: [v.copy(), v.copy()] for k, v in self.variables.items()} # add two in case of the first step

        
        self.lengths = [
            len(self.variables['v_p']),
            len(self.variables['v_a']),
            len(self.variables['lambda']),
            len(self.variables['q']),
        ]
        
        for k, v in self.variables.items():
            print(k, v.shape, sep='\t')

        
        all_blocks = [[cp.sparse.csr_matrix((self.lengths[i], self.lengths[j])) for j in range(4)] for i in range(4)]

        all_blocks[0][0] = cp.sparse.identity(len(self.variables['v_p']))

        all_blocks[1][1] = self.constants['M_ah'] / self.constants['delta_t']

        all_blocks[2][0] = 2 * self.constants['B_w_h']

        all_blocks[2][1] = 2 * -self.constants['B_gamma_h']
        all_blocks[1][2] = 2 * -constants['B_gamma_h'].T

        all_blocks[3][0] = -constants['J_h'].T
        all_blocks[3][3] = constants['M_qh'].T


        self.constants[r'\cos\sqrt K_h\Delta t'] = 0
        self.constants[r'\sin\sqrt K_h\Delta t/\sqrt K_h'] = 0
        self.constants[r'\frac{1-\cos(\sqrt K_h \Delta t)}{\sqrt K_h}'] = 0
        for omega in self.constants['omega_n'][:self.constants['n_specturm']]:

            alpha = (self.constants['R_p'] + self.constants['eta_p'] * omega) / 2
            k = omega ** 2 - alpha ** 2
            theta = np.sqrt(k) * self.constants['delta_t']

            self.constants[r'\cos\sqrt K_h\Delta t'] += np.exp(-alpha * self.constants['delta_t']) * np.cos(theta)
            self.constants[r'\sin\sqrt K_h\Delta t/\sqrt K_h'] += np.exp(-alpha * self.constants['delta_t']) * np.sin(theta) / np.sqrt(k)
            self.constants[r'\frac{1-\cos(\sqrt K_h \Delta t)}{\sqrt K_h}'] += np.exp(-alpha * self.constants['delta_t']) * (1 - np.cos(theta)) / np.sqrt(k)


        self.all_blocks = cp.sparse.bmat(all_blocks)
        assert self.all_blocks.shape == (sum(self.lengths), sum(self.lengths))
        # self.all_blocks = cupyx.scipy.sparse.linalg.factorized(self.all_blocks)

        self.currentTime = 0
        

    def update(self):
        
        self.string.get_f_sh(self.currentTime)
        f_sh = cp.asarray(self.string.f_sh).reshape(-1, 1)
        # update v_s directly
        self.variables['v_s'] = self.variables['v_s'] + self.constants['M_sh_inverse'] @ (f_sh - self.constants['D_h'] @ self.variables['q']) * self.constants['delta_t']

        # update p directly
        self.variables['p_a'] = self.variables['p_a'] - self.constants['M_pah_inverse'] * self.constants['Gh'].T * self.variables['v_a'] * self.constants['delta_t']

        

        # solve for the others
        ## make the right coefficient vector
        # sp.vstack([
        #     sp.csc_matrix((20400, 1)),
        #     sp.csc_matrix((382500, 1)),
        #     sp.csc_matrix((44761, 1)),
        #     sp.csc_matrix((100, 1)),
        # ])
        b = np.vstack([
            self.constants[r'\cos\sqrt K_h\Delta t'] * self.variables['v_p'] + self.constants[r'\sin\sqrt K_h\Delta t/\sqrt K_h'] * self.variables['a_p'] + self.constants[r'\frac{1-\cos(\sqrt K_h \Delta t)}{\sqrt K_h}'] * (
                self.constants['B_w_h'].T / self.constants['delta_t'] @ self.variables['lambda'] + self.constants['J_h'] / self.constants['delta_t'] / 2 * self.history['q'][-2]),
            self.constants['M_ah'] / self.constants['delta_t'] @ self.variables['v_a'] + self.constants['Gh'] @ self.variables['p_a'],
            2 * self.constants['B_w_h'] @ self.variables['v_p'] - self.constants['B_gamma_h'] @ self.history['v_a'][-2],
            self.constants['M_qh'] @ self.history['q'][-1] + self.constants['D_h'].T @ self.variables['v_s']
        ])

        # x = self.all_blocks(b)
        x = cupyx.scipy.sparse.linalg.spsolve(self.all_blocks, b)
    
        # extract v_p, v_a, lambda, q
        self.variables['v_p'] = x[:self.lengths[0]]
        self.variables['v_a'] = x[self.lengths[0]:self.lengths[0] + self.lengths[1]]
        self.variables['lambda'] = x[self.lengths[0] + self.lengths[1]:self.lengths[0] + self.lengths[1] + self.lengths[2]]
        self.variables['q'] = x[self.lengths[0] + self.lengths[1] + self.lengths[2]:]

        # update histroy
        for k, v in self.variables.items():
            if not k.startswith('last'):
                self.history[k].append(v.copy())


        self.currentTime += self.constants['delta_t']


    def get_history(self):
        return self.history


if __name__ == '__main__':
    from models import String, SoundBoard, Air, Gamma
    import pickle
    
    constants: dict[str, sp.spmatrix | cp.sparse.spmatrix | np.ndarray | cp.ndarray]
    string: String
    soundboard: SoundBoard
    air: Air
    gamma: Gamma
    constants, (string, soundboard, air, gamma) = pickle.load(open('constants.pkl', 'rb'))


    simulation = PhysicalSimulation(constants, string, soundboard, air, gamma)

    # Update the simulation for 10 time steps
    for _ in range(10):
        simulation.update()

    history = simulation.get_history()

    # Print history shapes to verify
    for k, v in history.items():
        print(f"{k}: {len(v)} time steps, shape {v[0].shape}")
