#!/usr/bin/env python
"""
The gb module provide general banded matrix functionality.
"""

import numpy as np

def gb2ge(*args):
    """Converts LAPACK-like general banded matrices to dense storage.
    The general band matrix gb must be stored following LAPACK conventions.
    m is the number of rows, n the number of columns, and kl and ku
    are the number of sub- and super-diagonals, respectively.

    Invoke as any of
             ge = gb2ge (gb, m, n, kl, ku, ld)
             ge = gb2ge (gb, m, n, kl, ku)
             ge = gb2ge (gb, n, kl, ku)
             ge = gb2ge (gb, n, k)
    where the first form is for general rectangular matrices, the second
    for general square matrices, and the last for square matrices with
    the same number of sub- and super-diagonals.

    >>> m, n, kl, ku = 5, 5, 2, 1
    >>> gb = [   0., 11.,  21.,  31.,
    ...         12., 22.,  32.,  42.,
    ...         23., 33.,  43.,  53.,
    ...         34., 44.,  54.,   0.,
    ...         45., 55.,   0.,   0.]
    >>> gb2ge(gb, m, n, kl, ku, kl + 1 + ku)
    matrix([[ 11.,  12.,   0.,   0.,   0.],
            [ 21.,  22.,  23.,   0.,   0.],
            [ 31.,  32.,  33.,  34.,   0.],
            [  0.,  42.,  43.,  44.,  45.],
            [  0.,   0.,  53.,  54.,  55.]])
    >>> (   gb2ge(gb, m, n, kl, ku, kl + 1 + ku)
    ...  == gb2ge(gb, m, n, kl, ku             )).all()
    True
    >>> (   gb2ge(gb, m, n, kl, ku, kl + 1 + ku)
    ...  == gb2ge(gb,    n, kl, ku             )).all()
    True
    """

    nargin = len(args)
    if nargin == 6:
        gb, m, n, kl, ku, ld = args
    elif nargin == 5:
        gb, m, n, kl, ku = args
        ld = kl + 1 + ku;
    elif nargin == 4:
        gb, n, kl, ku = args
        m  = n
        ld = kl + 1 + ku;
    elif nargin == 3:
        gb, n, kl = args
        m  = n
        ku = kl;
        ld = kl + 1 + ku;
    else:
        raise ValueError("Invalid call to gb2ge: " + args)

    gb = np.array(gb).reshape(-1)
    ge = np.asmatrix(np.zeros((m, n)))
    for j in range(0, m):
        for i in range(max(0, j - ku), min(m, j + kl + 1)):
            ge[i,j] = gb[j*ld + ku + i - j];
    return ge


if __name__ == "__main__":
    import doctest
    doctest.testmod()
