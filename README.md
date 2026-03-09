# pyfesm 

Python tools for parsing and post-processing **OpenMX** electronic structure calculations. 

The project provides: 
- A structured parser for OpenMX output files 
- A clean separation between file I/O and physics post-processing 
- Support for SCF calculations, band structures, NEGF transport, Mulliken populations, and band unfolding 

The goal is to provide a lightweight and transparent interface for extracting numerical data from OpenMX calculations and performing custom post-processing in Python. 

--- 
## Philosophy 

### Parser module 
- Reads OpenMX output files 
- Stores structured numerical data 
- Performs unit normalization (Hartree → eV, Bohr → Å) 

### Physics module 
- Performs post-processing on parsed data 
- Implements analysis routines (DOS, transport, unfolding, etc.) 

--- 
## Features 

- **SCF and system data** 
- **Band structure mapping** 
- **Mulliken population analysis** 
- **NEGF transport support** 
- **Band unfolding data extraction** 

--- 
## Installation 

Clone the repository: ```bash git clone https://github.com/samlipton/pyfesm.git cd pyfesm ``` 

Install dependencies: ```bash pip install numpy scipy ``` 

--- 
## Basic Usage 

### OpenMX output file parser
```python from pyfesm.openmx.parser import OpenMX <br> calc = OpenMX("Si", path=".") <br> print(calc.Utot) # total energy (eV) <br> print(calc.Nk) # k-grid ```

### Band Structure map
```python (kx, ky, kz), Ek = calc.eigenvalues ``` 
Returns: 
- Unique k-grid axes 
- Eigenvalue array with shape `(nkx, nky, nkz, nbands)` 

### Density of States 
```python E, DOS = calc.DoS() ``` 

### NEGF Transmission 
```python E, T = calc.G0() ``` 
Returns: 
- Energy grid 
- Transmission array reshaped onto the transport k-grid 

---
## Units 

Internally normalized to: 
- **Energy:** eV 
- **Length:** angstroem 

--- 
## Code Structure 

``` pyfesm/ 
│ 
├── openmx/ 
│ ├── parser.py 
│ └── utils.py 
│ 
├── physics/ 
│ └── README.md ```


