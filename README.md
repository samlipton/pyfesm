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
- Performs unit normalization on:
  - Electronic structure
  - Mulliken population
  - NEGF transmissions

### Physics module 
- Performs post-processing on parsed data:
  - Energy filtering
  - Band unfolding interpolation

--- 
## Installation 

Clone the repository: ```bash git clone https://github.com/samlipton/pyfesm.git cd pyfesm ``` 

Install dependencies: ```bash pip install numpy scipy ``` 

--- 
## Basic Usage 

### Output file parser
```python from pyfesm.openmx.parser import OpenMX``` <br>
```python calc = OpenMX("Si", path=".")``` <br>
```python print(calc.Utot) # total energy (eV)``` <br>

### Electronic structure
```python (kx, ky, kz), Ek = calc.eigenvalues ``` <br>
Returns: 
- Unique k-grid axes 
- Eigenvalue array with shape `(nkx, nky, nkz, nbands)` 

### Density of states 
```python E, DOS = calc.DoS() ``` <br>
Returns: 
- Energy grid
- Density of states 

### NEGF transmission 
```python E, T = calc.G0() ``` <br>
Returns: 
- Energy grid 
- Transmission array 

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
│ └── physics.py
│ └── utils.py  
│
└── README.md ```


