#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
import functools
import hicstuff as hcs
from matplotlib import pyplot as plt
from scipy import sparse

load_raw_matrix = functools.partial(np.genfromtxt,
                                    skip_header=True,
                                    dtype=np.float64)


def raw_cols_to_sparse(M):
    n = np.amax(M[:, :-1]) + 1

    row = M[:, 0]
    col = M[:, 1]
    data = M[:, 2]
    S = sparse.coo_matrix((data, (row, col)), shape=(n, n))
    return S


def sparse_to_dense(M):

    D = M.todense()
    E = D + np.transpose(D) - 2 * np.diag(np.diag(D))
    return E


def plot_matrix(M, vmax=99):
    plt.figure()
    plt.imshow(M, vmax=np.percentile(M, vmax),
               cmap='Reds', interpolation='none')
    plt.colorbar()
    plt.axis('off')
    plt.show()


if __name__ == "__main__":
    try:
        filename = sys.argv[1]
    except IndexError:
        print("Usage: ./vizmap.py abs_fragments_contacts_weighted.txt")
        quit()

    plot_matrix(sparse_to_dense(raw_cols_to_sparse(load_raw_matrix(filename))))