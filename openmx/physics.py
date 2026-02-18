"""
Physics structure utilities.

Pure numerical analysis layer.
No file parsing here.
"""

import numpy as np
from .parser import OpenMX
from scipy.constants import value


# =============================================================================
# Physical constants
# =============================================================================

a0 = value("Bohr radius") * 1e10
Ha = value("Hartree energy in eV")


# =============================================================================
# Energy window filtering
# =============================================================================

def filter_energy_window(data: OpenMX, Emin, Emax):
    """
    Return masked eigenvalues within energy window.
    """

    _,energies = data.eigenvalues()

    mask = (energies >= Emin) & (energies <= Emax)
    return np.where(mask, data.eigenvalues, np.nan)


# =============================================================================
# Fatbands (Lorentzian spectral weight)
# =============================================================================

def fatbands(k_orb, E_orb, proj, dE=0.01, DE=0.04, Dk=1e-2):
    """
    Compute spectral weight intensity map.
    """

    k = np.unique(k_orb)
    Emin, Emax = E_orb.min(), E_orb.max()
    E = np.arange(Emin, Emax + dE, dE)

    nk, nE = len(k), len(E)
    norb = proj.shape[1]
    Imap = np.zeros((nk, nE, norb))

    for i in range(len(k_orb)):
        dk = k[:, None] - k_orb[i]
        dE_ = E[None, :] - E_orb[i]

        L = 1.0 / ((dk / Dk) ** 2 + (dE_ / DE) ** 2 + 1)

        for j in range(norb):
            Imap[..., j] += L * proj[i, j]

    return k, E, Imap