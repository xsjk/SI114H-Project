import cupy as cp
import scipy.sparse as sp
import numpy as np
import cupyx.scipy.sparse.linalg
from tqdm import tqdm


def gpu_to_cpu(constants):
    for k, v in constants.items():
        if isinstance(v, cp.sparse.spmatrix):
            constants[k] = v.get()
        elif isinstance(v, cp.ndarray):
            constants[k] = cp.asnumpy(v)
    cp.get_default_memory_pool().free_all_blocks()


def cpu_to_gpu(constants):
    for k, v in constants.items():
        if isinstance(v, sp.spmatrix):
            constants[k] = cp.sparse.csr_matrix(v)
        elif isinstance(v, np.ndarray):
            constants[k] = cp.asarray(v)


def inv(A: cp.sparse.spmatrix, batch_size: int = 1000) -> cp.sparse.spmatrix:
    assert A.shape[0] == A.shape[1]
    assert batch_size > 0

    n = A.shape[0]

    columns = cp.empty((0,), dtype=cp.int32)
    rows = cp.empty((0,), dtype=cp.int32)
    values = cp.empty((0,), dtype=A.dtype)

    A_factorized = cupyx.scipy.sparse.linalg.factorized(A)
    for i in tqdm(range(0, n, batch_size)):
        b = cp.zeros((n, batch_size))
        for j in range(min(batch_size, n-i)):
            b[i+j, j] = 1
        x = A_factorized(b)

        nonzeros = cp.asnumpy(x).nonzero()
        nonzeros_gpu = cp.asarray(nonzeros)
        columns = cp.concatenate((columns, nonzeros_gpu[0]), axis=0)
        rows = cp.concatenate((rows, nonzeros_gpu[1] + i), axis=0)
        values = cp.concatenate((values, x[nonzeros]), axis=0)

    return cupyx.scipy.sparse.coo_matrix((values, (rows, columns)), shape=(n, n))


def set_tol(A: cp.sparse.spmatrix | sp.spmatrix, tol: float) -> None:
    if isinstance(A, sp.coo_matrix):
        A.data[np.abs(A.data) < tol] = 0
    elif isinstance(A, cp.sparse.coo_matrix):
        A.data[cp.abs(A.data) < tol] = 0
    else:
        raise ValueError(f"Invalid matrix type: {type(A)}")
    A.eliminate_zeros()