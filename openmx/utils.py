"""
Utilities routines.

Responsible only handling IO operations.
"""

import numpy as np
from pathlib import Path
from typing import Tuple

# =============================================================================
# Utility functions
# =============================================================================

def read_column(path: Path, index: int = 0, dtype=float) -> np.ndarray:
    data = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) > index:
                data.append(parts[index])
    return np.asarray(data, dtype=dtype)


def read_two_columns(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    E, V = [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2 and ['#'] not in parts[0]:
                E.append(parts[0])
                V.append(parts[1])
    return np.asarray(E, float), np.asarray(V, float)