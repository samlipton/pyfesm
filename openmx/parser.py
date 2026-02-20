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
    MULLIKEN = auto()
    FORCES = auto()    

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
    basis: dict = field(default_factory=dict)
    Na: int = 0
    atm: List[str] = field(default_factory=list)
    acell: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    gcell: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    ngrid: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=int))

    # SCF
    solver: str = ""
    Nk: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=int))
    T: float = 0.0

    ## NEGF 
    Nk_negf: np.ndarray = field(default_factory=lambda: np.zeros(2, dtype=int))
    Nk_tran: np.ndarray = field(default_factory=lambda: np.zeros(2, dtype=int))
    NE_tran: int = 0
    En_tran: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))
    _Eb_tran: Optional[np.ndarray] = None
    _eta: float = 1e-3
    mu_ll: float = 0.0
    mu_rl: float = 0.0

    # Energies
    Utot: float = 0.0
    Ekin: float = 0.0
    mu: float = 0.0
    
    ## Orbitals unfolding
    Eb_orbs: np.ndarray = field(default_factory=lambda: np.zeros(2))
    Nkt_orb: int = 0
    kp_orb: dict = field(default_factory=dict) # kpoints in path definition
    k_orb: np.ndarray = field(default_factory=lambda: np.zeros((0,3)))
    orb_id: List[List[int]] = field(default_factory=list)

    # Bands
    nbands: int = 0
    _kpoints: Optional[np.ndarray] = None
    _eigenvalues: Optional[List[List[float]]] = None

    # Geometry
    coord: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    force: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))

    # Mulliken
    MC: np.ndarray = field(default_factory=lambda: np.zeros((0, 2)))
    MC_orb: List[List[float]] = field(default_factory=list)

    def __post_init__(self):
        self.path = Path(self.path)
        if self.auto_read:
            self.read()

        # post-computations
        if self.NE_tran>0:
            self.En_tran = np.linspace(*self._Eb_tran,self.NE_tran)

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

                # to implement: basis set 

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
                # NEGF
                # ----------------------------------------------------------------

                elif key == "NEGF.scf.Kgrid":
                    self.Nk_negf = np.array(parts[1:3], dtype=int)

                elif key == "NEGF.tran.Kgrid":
                    self.Nk_tran = np.array(parts[1:3], dtype=int)

                elif key == "NEGF.tran.energydiv":
                    self.NE_tran = int(parts[1])

                elif key == "NEGF.tran.energyrange":
                    self._Eb_tran = np.array(parts[1:3], float)
                    self._eta = float(parts[3])

                elif parts[:5] == ["Chemical","potential","of","left","lead"]:
                    self.mu_ll = float(parts[6]) * Ha

                elif parts[:5] == ["Chemical","potential","of","right","lead"]:
                    self.mu_rl = float(parts[6]) * Ha

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
                    nk = np.prod(self.Nk)
                    self._kpoints = np.zeros((nk, 3))
                    self._eigenvalues = []
                    self._ik = 0

                elif parts[::2] == ["k1=", "k2=", "k3="]:
                    state = ParseState.EIGENVALUES
                    self._kpoints[self._ik] = np.array(parts[1::2], float)
                    self._eigenvalues.append([])

                elif state == ParseState.EIGENVALUES:
                    if not key.isdigit(): # or band_index < nbands_max
                        state = ParseState.NONE
                        continue 

                    # band_index = int(parts[0]) - 1
                    energy = float(parts[1]) * Ha
                    self._eigenvalues[self._ik].append(energy)

                # -------------------------------------------------------------
                # Mulliken population
                # -------------------------------------------------------------

                elif parts in (
                    ["Up", "spin", "Down", "spin", "Sum", "Diff"],
                    ["Up", "Down", "Sum", "Diff", "theta", "phi"],
                    ):
                    state = ParseState.MULLIKEN
                    self.MC = np.zeros((self.Na, 2))

                    for i in range(self.Na):
                        parts = next(f).split()
                        self.MC[i] = list(map(float, parts[2:4]))

                    state = ParseState.NONE

                elif parts == ["Decomposed", "Mulliken", "populations"]:
                    for atom in self.atm:
                        basis_str = self.basis[atom]

                        # Extract orbital counts from basis string
                        # e.g. "C5.0-s2p2d1"
                        counts = basis_str.split("-")[1]
                        s, p, d = 0, 0, 0
                        if "s" in counts:
                            s = int(counts.split("s")[1][0])
                        if "p" in counts:
                            p = int(counts.split("p")[1][0])
                        if "d" in counts:
                            d = int(counts.split("d")[1][0])

                        s_list, p_list, d_list = [], [], []
                        nlines = 2*s + 4*p + 6*d + 6
                        for _ in range(nlines):
                            parts = next(f).split()
                            if not parts:
                                continue
                            if parts[0] == "s":
                                s_list.append(parts[2:4])
                            elif parts[0] in ("px","py","pz"):
                                p_list.append(parts[2:4])
                            elif parts[0] in ("d3z^2-r^2","dx^2-y^2","dxy","dxz","dyz"):
                                d_list.append(parts[2:4])
                        t_list = s_list + p_list + d_list 
                        self.MC_orb.append(t_list)

                # -------------------------------------------------------------
                # Bands unfolding
                # -------------------------------------------------------------

                elif key == "Unfolding.LowerBound":
                    self.Ebounds[0] = float(parts[1])

                elif key == "Unfolding.UpperBound":
                    self.Ebounds[1] = float(parts[1])

                elif key == "Unfolding.desired_totalnkpt":
                    self.Nkpoints = int(parts[1])

                elif parts == ["Unfolding", "calculation", "for", "band", "structure"]:
                    self._parse_unfolding_block(f)                

                elif parts[:3] == ["For", "each", "path:"]: # to implement
                    vals = list(map(int, parts[4:4+self.Nkpoints+1]))
                    #self.Nkpath_orb = np.array(vals)
                    #self.ikpath_orb = np.cumsum([0] + vals[:-1])

                # to implement: ka kb kc and k_orb

                # -------------------------------------------------------------
                # Forces
                # -------------------------------------------------------------

                elif key == "<coordinates.forces":
                    state = ParseState.FORCES
                    next(f)

                    self.coord = np.zeros((self.Na, 3))
                    self.force = np.zeros((self.Na, 3))

                    for i in range(self.Na):
                        line = next(f).split()
                        vals = list(map(float, line[2:]))
                        self.coord[i] = vals[:3]
                        self.force[i] = vals[3:]

                    state = ParseState.NONE

    def _parse_unfolding_block(self, f): # to implement
        # skip header
        for _ in range(11):
            next(f)

        return None

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

    def G0(self, grid=None):

        if grid is None:
            grid = self.Nk[:2]

        T = []
        E = None

        for i in range(grid[0]):
            for j in range(grid[1]):
                file = self.path / f"{self.name}.tran{i}_{j}"
                Ei, Ti = read_two_columns(file)
                if E is None:
                    E = Ei
                T.append(Ti)

        T = np.array(T).reshape(grid[0], grid[1], -1)

        return E, T
