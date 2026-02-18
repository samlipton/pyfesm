"""
Data model and parser for OpenMX calculations.

returns structured OpenMX data after reading output files
"""

from __future__ import annotations
import numpy as np
from .utils import read_two_columns
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import Tuple, List, Optional
from scipy.constants import value


# =============================================================================
# Physical constants
# =============================================================================

a0 = value("Bohr radius") * 1e10
Ha = value("Hartree energy in eV")


# =============================================================================
# Parser state machine
# =============================================================================

class ParseState(Enum):
    NONE = auto()
    EIGENVALUES = auto()


# =============================================================================
# OpenMX Data Container
# =============================================================================

@dataclass
class OpenMX:

    name: str
    path: Path = Path(".")
    ext: str = "out"
    auto_read: bool = True

    # System
    Ns: int = 0
    species: List[str] = field(default_factory=list)
    Na: int = 0
    atm: List[str] = field(default_factory=list)
    acell: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    gcell: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    ngrid: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=int))

    # SCF
    solver: str = ""
    Nk: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=int))
    T: float = 0.0

    # Energies
    Utot: float = 0.0
    Ekin: float = 0.0
    mu: float = 0.0

    # Bands
    # nbands: int = 0
    _kpoints: Optional[np.ndarray] = None
    _eigenvalues: Optional[np.ndarray] = None

    # Geometry
    coord: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    force: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))

    # Mulliken
    MC: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))

    def __post_init__(self):
        self.path = Path(self.path)
        if self.auto_read:
            self.read()

    # =========================================================================
    # Main parser
    # =========================================================================

    def read(self) -> None:

        state = ParseState.NONE
        file_path = self.path / f"{self.name}.{self.ext}"

        with open(file_path) as f:
            for raw in f:
                parts = raw.split()
                if not parts:
                    continue

                key = parts[0]

                # ----------------------------------------------------------------
                # System block
                # ----------------------------------------------------------------

                if key == "Species.Number":
                    self.Ns = int(parts[1])

                elif key == "Atoms.Number":
                    self.Na = int(parts[1])

                elif key == "<Atoms.UnitVectors":
                    self.acell = np.array(
                        [np.fromstring(next(f), sep=" ") for _ in range(3)]
                    )

                elif parts[:3] == ["Cell", "vectors", "(bohr)"]:
                    self.gcell = np.array(
                        [np.fromstring(next(f)[12:], sep=",") for _ in range(3)]
                    )

                # ----------------------------------------------------------------
                # SCF
                # ----------------------------------------------------------------

                elif key == "scf.EigenvalueSolver":
                    self.solver = parts[1]

                elif key == "scf.Kgrid":
                    self.Nk = np.array(parts[1:4], dtype=int)

                elif key == "scf.ElectronicTemperature":
                    self.T = float(parts[1])

                # ----------------------------------------------------------------
                # Energies
                # ----------------------------------------------------------------

                elif key == "Utot.":
                    self.Utot = float(parts[1]) * Ha

                elif key == "Ukin.":
                    self.Ekin = float(parts[1]) * Ha

                elif parts[:2] == ["Chemical", "Potential"]:
                    self.mu = float(parts[4]) * Ha

                # ----------------------------------------------------------------
                # Band section detection
                # ----------------------------------------------------------------

                elif parts == ["Eigenvalues","(Hartree)","of","SCF","KS-eq."]:
                    self._init_band_storage()

                elif parts[::2] == ["k1=", "k2=", "k3="]:
                    state = ParseState.EIGENVALUES
                    self._parse_kpoint(parts)

                elif state == ParseState.EIGENVALUES:
                    if self._parse_eigen_line(parts):
                        state = ParseState.NONE

    # =========================================================================
    # Band parsing internals
    # =========================================================================

    def _init_band_storage(self):
        nk = np.prod(self.Nk)
        self._kpoints = np.zeros((nk, 3))
        self._eigenvalues = []
        self._ik = 0

    def _parse_kpoint(self, parts):
        self._kpoints[self._ik] = np.array(parts[1::2], float)
        self._eigenvalues.append([])

    def _parse_eigen_line(self, parts):
        if not parts[0].isdigit(): # or band_index < nbands_max
            self._ik += 1
            return True

        # band_index = int(parts[0]) - 1
        energy = float(parts[1]) * Ha

        self._eigenvalues[self._ik].append(energy)

    # =========================================================================
    # Public band interface
    # =========================================================================

    @property
    def eigenvalues(self) -> Tuple[np.ndarray, np.ndarray]:

        if self._eigenvalues is None:
            raise RuntimeError("Eigenvalues not parsed.")

        k = self._kpoints
        Ek = np.array(self._eigenvalues)

        # reshape onto 3D grid
        kx, ix = np.unique(k[:, 0], return_inverse=True)
        ky, iy = np.unique(k[:, 1], return_inverse=True)
        kz, iz = np.unique(k[:, 2], return_inverse=True)

        return (kx, ky, kz), Ek

    # =========================================================================
    # DOS
    # =========================================================================

    def DoS(self, subdir="PDoS") -> Tuple[np.ndarray, np.ndarray]:
        solver_name = "NEGF" if self.solver == "NEGF" else "Tetrahedron"
        file = self.path / subdir / f"{self.name}.DOS.{solver_name}"
        return read_two_columns(file)
    
    # =========================================================================
    # Transmission
    # =========================================================================

    def G0(self, grid=None, path=None):

        if grid is None:
            grid = self.Nk[:2]

        if path is None:
            path = self.path

        T = []
        E = None

        for i in range(grid[0]):
            for j in range(grid[1]):
                file = path / f"{self.name}.tran{i}_{j}"
                Ei, Ti = read_two_columns(file)
                if E is None:
                    E = Ei
                T.append(Ti)

        T = np.array(T).reshape(grid[0], grid[1], -1)

        return E, T