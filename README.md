# GRBoondi

GRBoondi is an open-sourced code for simulating generalized Proca fields on arbitrary analytic fixed backgrounds, based on the publicly available 3+1D numerical relativity code GRChombo. The goal of GRBoondi is to reduce the prequisite knowledge of numerical relativity and GRChombo in the numerical studies of generalized Proca theories . The main steps for performing a study should only be inputting the additions to the equations of motion beyond the base Proca theory and GRBoondi can automatically incorporate the higher-order terms in the simulation.



GRBoondi is written entirely in C++14, using hybrid MPI/OpenMP 
parallelism to achieve good performance on the latest architectures.
It inherits all of the capabilities of the main [GRChombo](https://github.com/GRChombo/GRChombo) 
code, which makes use of the Chombo library for adaptive mesh refinement.

## Other naming ideas

[Boondi](https://drive.google.com/file/d/1fTf8XjkJ62fGTQJdoImrA8178-2kyxeG/view): Koorie word for Multipurpose tool weapon, hunting implement

[Milijun](https://aphasialab.org/wagiman/dict/dict.html): Wagiman word for morning star 

[Pirri](https://australianwords.au/advanced-search/): Arabana word for tool

[Woomera](https://en.wikipedia.org/wiki/Woomera_(spear-thrower)): Aboriginal wooden throwing spear

[Boonta](https://drive.google.com/file/d/1fTf8XjkJ62fGTQJdoImrA8178-2kyxeG/view) Mentally unwell, crazy (from the Koorie tribe)

## Getting started
WIP

## Contributing
We welcome feedback, bug reports, and contributions. Please consult the repo owner for contributions and feel free to open pull requests. For coding style, please see the GRChombo repo wiki.

## License
GRBoondi is licensed under the BSD 3-Clause License. Please see LICENSE for details.

## Citation
WIP