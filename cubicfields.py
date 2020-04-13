import numpy as np
from functools import lru_cache
from numpy.linalg import matrix_power as matpow


def jz(j):
    """ Returns the jz matrix for a given angular momentum j."""
    return np.diag(np.arange(-j, j + 1))


def _jpcoeff(l, m):
    """ Helper function for jplus and jminus."""
    return np.sqrt(l * (l + 1) - m * (m + 1))


def jp(j):
    """ Returns the jplus matrix for a given angular momentum j."""
    return np.diag(_jpcoeff(j, np.arange(-j, j)), k=1)


def jm(j):
    """ Returns the jminus matrix for a given angular momentum j."""
    return np.diag(_jpcoeff(j, np.arange(-j, j)), k=-1)


def ident(j):
    return np.eye(int(2 * j + 1))


def o_4_0(j):
    """ Returns the o_4_0 matrix for a given angular momentum j."""
    return (
        # fmt: off
        35 * matpow(jz(j), 4)
        - (30 * j * (j + 1) - 25) * matpow(jz(j), 2)
        - 6 * j * (j + 1) * ident(j)
        + 3 * j**2 * (j + 1)**2 * ident(j)
        # fmt: on
    )


def o_4_4(j):
    """ Returns the o_4_4 matrix for a given angular momentum j."""
    return 0.5 * (matpow(jp(j), 4) + matpow(jm(j), 4))


def o_6_0(j):
    """ Returns the o_6_0 matrix for a given angular momentum j."""

    return (
        # fmt: off
        231 * matpow(jz(j), 6)
        - 105 * (3 * j * (j + 1) - 7) * matpow(jz(j), 4)
        + (105 * j**2 * (j + 1)**2 - 525 * j * (j + 1) + 294) * matpow(jz(j), 2)
        - 5 * j**3 * (j + 1)**3 * ident(j)
        + 40 * j**2 * (j + 1)**2 * ident(j)
        - 60 * j * (j + 1) * ident(j)
        # fmt: on
    )


def o_6_4(j):
    """ Returns the o_6_6 matrix for a given angular momentum j."""
    t1 = 11 * matpow(jz(j), 2) - j * (j + 1) * ident(j) - 38 * ident(j)
    t2 = matpow(jp(j), 4) + matpow(jm(j), 4)
    return 0.25 * (t1 @ t2 + t2 @ t1)


@lru_cache(maxsize=4)
def get_o_4(j):
    """ Returns o_4 = o_4_0 + 5 * o_4_4."""
    return o_4_0(j) + 5 * o_4_4(j)


@lru_cache(maxsize=4)
def get_o_6(j):
    """ Returns o_6 = o_6_0 - 21 * o_6_4."""
    return o_6_0(j) - 21 * o_6_4(j)


def energies(j, W, x, f4=60, f6=13860):
    """ Returns a list of the eigenenergies of the Hamiltonian specified
    by the total angular momentum j, the scale parameter W, and the mixing
    parameter x.
    """
    if j == 13 / 2:
        f4 = 60
        f6 = 7560
    if j == 15 / 2:
        f4 = 60
        f6 = 13860
    return W * np.linalg.eigvalsh(x * get_o_4(j) / f4 + (1 - abs(x)) * get_o_6(j) / f6)
