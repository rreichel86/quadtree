
# A Quadtree based finite element mesh generator to discretize two-dimensional heterogeneous solids

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

The two-dimensional heterogeneous solid is given by its outer boundary
and material interfaces. These are later modeled by simple closed polygonal
chains.  

The Figures (a) to (j) illustrate the Quadtree subdivision of a circular disk
with a circular hole and a circular disk with a circular inclusion. The
corresponding finite element meshes are shown in Figures (k) and (l).

| | | |
| :---: | :---: | :---: |
| <img src="./Images/Homogen.png" alt=" " width="150px"/> | <img src="./Images/Heterogen.png" alt=" " width="150px"/>  | <img src="./Images/QtreeBSP_0.png" alt=" " width="150px"/> |
| (a)  | (b)  | (c) |
| <img src="./Images/QtreeBSP_1.png" alt=" " width="150px"/> | <img src="./Images/QtreeBSP_2.png" alt=" " width="150px"/> | <img src="./Images/QtreeBSP_3.png" alt=" " width="150px"/> |
| (d) | (e) | (f) |
| <img src="./Images/QtreeBSP_4.png" alt=" " width="150px"/> | <img src="./Images/QtreeBSP_5.png" alt=" " width="150px"/> | <img src="./Images/QtreeBSP_6.png" alt=" " width="150px"/> |
| (g) | (h) | (i) |
| <img src="./Images/QtreeBSP_7.png" alt=" " width="150px"/> | <img src="./Images/QtreeBSP_8.png" alt=" " width="150px"/> | <img src="./Images/QtreeBSP_9.png" alt=" " width="150px"/> |
| (j) | (k) | (l) |




### Table of contents  

1. [What's in the directory? ](#whats-in-the-directory)
2. [Using the Quadtree mesh generator](#using-the-quadtree-mesh-generator)
3. [Acknowledgments](#acknowledgments)
4. [Contact](#contact)

## What's in the directory? <a name="whats-in-the-directory"></a>

.  
├── Documentation  
├── Examples  
├── Makefile  
├── Output/mesh  
├── Output/plot  
├── QtreePlotMesh.m  
├── README.md  
└── src  

| File            | Description |
| :-------------- | :---- |
| [Documentation](./Documentation/) | Documentation |
| [Examples](./Examples/)            | Contains examples of how to define the domain to be meshed |
| Makefile            | Makefile to build the Quadtree mesh generator |
| Output/mesh         | polygonal mesh data |
|                     | `scor.txt`  - File containing the nodal coordiantes  |
|                     | `selm.txt` - File containing the element connectivity |
| Output/plot         | polygonal mesh plot as a PNG image |
|                     | `QuadtreeMesh.png` (default name)|
| QtreePlotMesh.m     | MATLAB script to plot the Quadtree mesh |
| README.md           | This file |
| src                 | Source code of the Quadtree mesh generator|


## Using the Quadtree mesh generator

### Command Line Options

```
 -h or --help           Print this help message and exit
 -i or --iteractive     Iteractive mode

```

To run, for example, the [YetiFootprint.txt](./Examples/YetiFootprint.txt) example in Examples/
```
 $ ./Quadtree Examples/YetiFootprint.txt

```
The polygonal mesh data is stored in ./Output/mesh/scor.txt and
./Output/mesh/selm.txt; `scor.txt` contains the nodal coordinates and
`selm.txt` the element connectivity.


## Acknowledgments <a name="acknowledgments"></a>

The financial support of the DFG (German Research Foundation) under 
[project number 433734847](https://gepris.dfg.de/gepris/projekt/433734847?language=en)
is gratefully acknowledged.

## Contact <a name="contact"></a>  


**Rainer Reichel** M.Sc.  
LBB - Lehrstuhl für Baustatik und Baudynamik  
RWTH Aachen University 
email: <reichel@lbb.rwth-aaachen.de>


