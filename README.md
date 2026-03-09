# pyfesm 

PYthon tools For Electronic Structure Calculations. 

The project provides: 
- Parsing for output files 
- Post-processing data for physical analysis 

The goal is to provide a lightweight and transparent interface in python for the **OpenMX** software package. 

--- 
## Implementation

### Parser module 
- Reads OpenMX output files 
- Stores structured numerical data 
- Performs unit normalization (Hartree → eV, Bohr → Å) 

### Physics module 
- Performs post-processing on parsed data 
 
#### Features 

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
```python from pyfesm.openmx.parser import OpenMX```
```python calc = OpenMX("Si", path=".")```
```python print(calc.Utot) # total energy (eV)```

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


