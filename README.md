# Research Code for _Quantum State Transfer in Interacting, Multiple-Excitation Systems_

This repository contains the code used in the research paper "Quantum State Transfer in Interacting, Multiple-Excitation Systems" https://arxiv.org/abs/2405.06853. This repo contains all the codes, primarily written in Julia and Python, that were used for data analysis, simulations, and optimizations discussed in the paper.

### Overview

This repository includes several scripts and notebooks used during the research for "Quantum State Transfer in Interacting, Multiple-Excitation Systems." The key feature of this work is the optimization algorithm, "Dual Annealing" used to locate sets of couplings that permit near perfect QST in systems with multiple interacting excitaitons. 

The repository is organized as follows:
```
multi_exc_QST_notebook.ipynb 
- A Julia Jupyter Notebook containing the Dual Annealing code, the basis and Hamiltonian creation, the time evolution simulation, and various analysis elements.

PLTQSTGraph.ipynb
- A Python Jupyter Notebook used to create all of the plots and figures featured in the paper.

QSTData.jl
- The Julia script file that was run on a computer cluster to benchmark the convergence of our optimization methods.
```

### Installation

To run the code provided in this repository, you need to have Python 3 with the ipykernal installed and Julia (1.7.3) with the iJulia kernal installed. The necessary dependancies for each notebook are listed in the first cell. 

### License

This repository is licensed under the MIT License - see the LICENSE file for details.

### Contact

For any questions or inquiries regarding the code or the research, please contact Alexander Yue at alexyue@stanford.edu.
