
# A Quadtree based finite element mesh generator to discretize two-dimensional heterogeneous solids

To discretize two-dimensional heterogeneous solids, a quadtree based finite
element mesh generator is developed. The solid is given by its outer boundary
and material interfaces. These are later modeled by simple closed polygonal
chains. 

| | | |
| :---: | :---: | :---: |
| <img src="./Images/Homogen.png" alt=" " width="150px"/> | <img src="./Images/Heterogen.png" alt=" " width="150px"/>  | <img src="./Images/QtreeBSP_0.png" alt=" " width="150px"/> |
| <img src="./Images/QtreeBSP_1.png" alt=" " width="150px"/> | <img src="./Images/QtreeBSP_2.png" alt=" " width="150px"/> | <img src="./Images/QtreeBSP_3.png" alt=" " width="150px"/> |
| <img src="./Images/QtreeBSP_4.png" alt=" " width="150px"/> | <img src="./Images/QtreeBSP_5.png" alt=" " width="150px"/> | <img src="./Images/QtreeBSP_6.png" alt=" " width="150px"/> |
| <img src="./Images/QtreeBSP_7.png" alt=" " width="150px"/> | <img src="./Images/QtreeBSP_8.png" alt=" " width="150px"/> | <img src="./Images/QtreeBSP_9.png" alt=" " width="150px"/> |


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
├── QtreePlotMesh.m  
├── README.md  
└── src  

| File            | Description |
| :-------------- | :---- |
| [Documentation](./Documentation/) | Documentation |
| [Examples](./Examples/)            | Contains examples of how to define the domain to be meshed |
| Makefile            | Makefile to build the Quadtree mesh generator |
| QtreePlotMesh.m     | MATLAB script to plot the Quadtree mesh |
| README.md           | This file |
| src                 | Source code of the Quadtree mesh generator|


## Using the Quadtree mesh generator

### Command Line Options

```
 -h or --help           Print this help message and exit
 -i or --iteractive     Iteractive mode

```

To run the [YetiFootprint.txt](./Examples/YetiFootprint.txt) example in Examples/
```
 $ ./Quadtree Examples/YetiFootprint.txt

```


## Acknowledgments <a name="acknowledgments"></a>

The financial support of the DFG (German Research Foundation) under 
[project number 433734847](https://gepris.dfg.de/gepris/projekt/433734847?language=en)
is gratefully acknowledged.

## Contact <a name="contact"></a>  


**Rainer Reichel** M.Sc.  
LBB - Lehrstuhl für Baustatik und Baudynamik  
RWTH Aachen University 
email: <reichel@lbb.rwth-aaachen.de>


