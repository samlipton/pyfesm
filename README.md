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
  - SCF and System data
  - Electronic band structure
  - Mulliken population
  - NEGF transmissions

### Physics module 
- Performs post-processing on parsed data:
  - Energy filtering
  - Band unfolding interpolation

--- 
## Installation 

Clone the repository: ```$git clone https://github.com/samlipton/pyfesm.git cd pyfesm ``` 

Install dependencies: ```$pip install numpy scipy ``` 

--- 
## Basic Usage 

### Output file parser
```from pyfesm.openmx.parser import OpenMX``` <br>
```calc = OpenMX("Si", path=".")``` <br>
```print(calc.Utot) # total energy (eV)``` <br>

### Electronic band structure
```python (kx, ky, kz), Ek = calc.eigenvalues ``` <br>
Returns the k-grid axes Eigenvalue array with shape `(nkx, nky, nkz, nbands)` 

### Density of states 
```python E, DOS = calc.DoS() ``` <br>
Returns the Energy grid and Density of states 

### NEGF transmission 
```python E, T = calc.G0() ``` <br>
Returns the Energy grid and Transmission array 

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


