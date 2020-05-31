from scipy.sparse.coo import coo_matrix
import globals as g
import numpy as np
import _fwi_ls


def _build_extended_A(A_row, A_col, A_data, n_rows, n_cols):
    big_A_coo_i = [None] * (4 * len(A_row))
    big_A_coo_j = [None] * (4 * len(A_row))
    big_A_coo_data = [None] * (4 * len(A_row))
    for idx, (i, j, data) in enumerate(zip(A_row, A_col, A_data)):
        big_A_coo_i[4 * idx] = i
        big_A_coo_j[4 * idx] = j
        big_A_coo_data[4 * idx] = data.real

        big_A_coo_i[4 * idx + 1] = i + n_rows
        big_A_coo_j[4 * idx + 1] = j
        big_A_coo_data[4 * idx + 1] = data.imag

        big_A_coo_i[4 * idx + 2] = i
        big_A_coo_j[4 * idx + 2] = j + n_cols
        big_A_coo_data[4 * idx + 2] = -data.imag

        big_A_coo_i[4 * idx + 3] = i + n_rows
        big_A_coo_j[4 * idx + 3] = j + n_cols
        big_A_coo_data[4 * idx + 3] = data.real
    return big_A_coo_i, big_A_coo_j, big_A_coo_data


def convert_complex_linear_system_to_real(A, b):
    build_extended_A = g.choose_impl(_fwi_ls.build_extended_A, _build_extended_A)

    big_A_coo_i, big_A_coo_j, big_A_coo_data = build_extended_A(
        A.row, A.col, A.data, A.shape[0], A.shape[1],
    )

    big_A = coo_matrix(
        (big_A_coo_data, (big_A_coo_i, big_A_coo_j)),
        shape=(2 * A.shape[0], 2 * A.shape[1]),
        dtype=np.float,
    )

    # b vector is not sparse
    big_b = np.zeros(shape=(2 * len(b), 1))
    big_b[0 : len(b)] = b.real
    big_b[len(b) : 2 * len(b)] = b.imag

    return big_A, big_b


def extract_complex_solution(big_x):
    x = np.zeros(shape=(len(big_x) // 2,), dtype=np.complex)
    x.real = big_x[0 : len(x)]
    x.imag = big_x[len(x) : 2 * len(x)]
    return x
