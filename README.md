# GRBoondi

[![status](https://joss.theoj.org/papers/9565689c5c5d5da3c39adb87a4e7d255/status.svg)](https://joss.theoj.org/papers/9565689c5c5d5da3c39adb87a4e7d255)

GRBoondi is an open-sourced code for simulating generalized Proca fields on arbitrary analytic fixed backgrounds, based on the publicly available 3+1D numerical relativity code GRChombo. The goal of GRBoondi is to reduce the prequisite knowledge of numerical relativity and GRChombo in the numerical studies of generalized Proca theories . The main steps for performing a study should only be inputting the additions to the equations of motion beyond the base Proca theory and GRBoondi can automatically incorporate the higher-order terms in the simulation.



GRBoondi is written entirely in C++14, using hybrid MPI/OpenMP 
parallelism to achieve good performance on the latest architectures.
It inherits all of the capabilities of the main [GRChombo](https://github.com/GRChombo/GRChombo) 
code, which makes use of the Chombo library for adaptive mesh refinement.


## Getting started
See the [wiki pages](https://github.com/ShaunFell/GRBoondi/wiki)

## Contributing
We welcome feedback, bug reports, and contributions. Please consult the repo owner for contributions and feel free to open pull requests. For coding style, please see the GRChombo repo wiki.

## License
GRBoondi is licensed under the BSD 3-Clause License. Please see LICENSE for details.

## Citation
This repository has been submitted to the Journal of Open Source Software. In the mean-time, if you use this software for your research, please cite our arXiv paper:

https://arxiv.org/abs/2405.01348