---
title: 'GRBoondi: A code for evolving Generalized Proca theories on arbitrary backgrounds'
tags:
  - c++
  - MPI
  - Open MP
  - vector intrinsics
  - Generalized Proca
  - gravity
  - general relativity
  - numerical relativity
  - adaptive refinement
authors:
- name: Shaun David Brocus Fell
  orcid: 0000-0002-8059-0359
  affiliation: 1
- name: Lavinia Heisenberg
  affiliation: 1
affiliations:
- name: Institute for Theoretical Physics, Universitaet Heidelberg, Philosophenweg 12, 69120 Heidelberg, Germany
  index: 1
date: 03-05-2024
bibliography: paper.bib
---

# Summary

Proca theories are a simple massive extension of the typical electromagnetic field. Such a simple modification of the underlying theory can have tremendous phenomenological implications, particularly in the case of early-universe high-energy physics, structure formation, and various dark matter models. On the theoretical side, they arise naturally in certain string theories and axionic interactions. Phenomenologically, Proca fields can be produced in the early universe from quantum fluctuations. This leads the Proca field to be an excellent candidate for the dark matter particle. Naturally, generalizations of the Proca field can lead to even richer astrophyical implications, providing possible solutions to many challenges across a wide range of astrophysical phenomena.

## Generalized Proca 

Generalized Proca theories [@Heisenberg_2019;@Heisenberg_2014] provide a rich landscape to search for solutions to deep, fundamental questions, such as the nature of dark energy and dark matter. Moreover, strong gravity regimes, like those surrounding dense, compact astrophysical entities, offer novel avenues to explore fundamental fields. To effectively probe these domains, precise models are imperative for sifting through vast quantities of data. Crafting concrete models for the entirety of generalized Proca theories represents a formidable endeavor, typically relying on numerical methods. These numerical methods are usually written to be case dependent, especially for a particular background. This makes it difficult to generalize the computational code to account for higher-order couplings or different types of backgrounds. The difficulty is especially amplified when considering the full landscape of generalized Proca theories and beyond [@Heisenberg:2016eld]. Already at the level of the Lagrangian, the equations are immense,
$$
	\mathcal{L}_{\text{g.P.}} = \sqrt{-g} \left(-\frac{1}{4}F_{\mu \nu}F^{\mu \nu} + \displaystyle\sum_{n=2}^{6} \alpha_n \mathcal{L}_n \right) \;,
$$
where each sub-Lagranian $\mathcal{L}_n$ contains terms of different orders in the derivatives of the Proca field [@Heisenberg_2019;@Heisenberg_2014]. Even including solely $\mathcal{L}_2$, the equations of motion can be extremely cumbersome. Solving these equations of motion analytically very quickly becomes intractable.

GRBoondi provides a unified interface for computing the evolution of any generalized Proca model on an arbitrary background. Given a specific background with known expressions for the metric variables and initial data for the Proca field, GRBoondi numerically computes the time evolution of the Proca field, given user-specified generalized Proca equations of motion. While the simultaneous evolution of the metric and matter fields offers the most comprehensive depiction of their temporal progression, there are instances where the density of the Proca field is negligible in comparison to the background curvature. In such cases, fixed background evolution routines serve as an excellent approximation to the complete mutual evolution. A core feature of GRBoondi is that it simplifies many of the boiler-plate code required to begin a simulation from scratch. The only additions a user needs to input are the initial conditions, the modifications to the equations of motion (EOM), and the background functions. GRBoondi will automatically compute various diagnostic quantities and plot files. Moreover, GRBoondi is incredibly modular and modifications to any of the in-built functions is effortless.

The fixed background approximation provides a huge speedup compared to the full evolution of the Einstein equations, since only a small fraction of the variables actually require a time integration. Not only do simulations now require far less resources, the speed up can be up to hundreds of times faster. GRBoondi is based on the publicly available numerical relativity (NR) code GRChombo [@Andrade2021;@Clough2015sqa], which itself is based on the open source adaptive mesh refinement (AMR)-based differential equation solver Chombo [@Adams:2015kgr]. Some smaller pieces of GRBoondi are also based off the recently released open source software GRDzhadzha [@Aurrekoetxea2024], which is also based off of GRChombo. These comprise the core dependencies of GRBoondi.

# Statement of need

In practice, any NR software library can be used to evolve generalized Proca theories. Prominent examples of NR codes encompass the extensive Einstein toolkit [@EinsteinToolkit:2023_11] and its related Cactus framework [@Goodale:2002a], Kranc [@Kranc:web], LEAN [@Sperhake:2006cy], and Canuda [@Canuda:zenodo]. Other libraries worth mentioning are the non-public BAM [@Bruegmann:2006ulg], AMSS-NCKU [@Galaviz:2010mx], PAMR [@East:2011aa] and HAD [@Neilsen:2007ua]. The non-exhaustive list can be expanded with SPeC [@Pfeiffer:2002wt], which is a pseudo-spectral code that uses generalized harmonic coordinates. Similarly, SpECTRE [@deppe_nils_2021_4734670;@Kidder:2016hev;@Cao:2018vhw] uses the Galerkin methods. NRPy [@Ruchlin:2017com] is a python library that aims for use on non-high performance computing clusters. There are also cosmological simulation codes like CosmoGRaPH [@Mertens:2015ttp] and GRAMSES [@Barrera-Hinojosa:2019mzo]. Simflowny [@Palenzuela:2018sly] is a magneto-hydrodynamic simulation software. GRAthena++ [@Daszuta:2021ecf] uses the oct-tree AMR for maximum scaling. ExaGRyPE is a numerical relativity program leveraging the ExaHyPRE PDE solver [@zhang2024exagrypenumericalgeneralrelativity]. GRDzhadzha [@Aurrekoetxea2024] is another fixed background code, based off the GRChombo [@Andrade2021] framework, which GRBoondi utilizes.

This extensive list underscores the abundance of NR libraries at one's disposal. However, none of them provide a tailored unified interface for studying the vast landscape of gravity theories, like generalized Proca. Using any of the numerous frameworks requires significant work in order to evolve even a single generalized Proca theory. GRBoondi tackles this problem by providing a collection of specialized tools for computing the generalized Proca equations, relying on existing tools which allows for rapid updating and debugging. On top of this, GRBoondi offers catered plotting routines for viewing data, leveraging the highly parallelizable VisIt [@visit_dav] analysis tool. In addition, since the metric variables and their derivatives are computed exactly at each grid point, the adaptability of the AMR grid can be focused solely on the matter variables. 

While backreaction is neglected in GRBoondi, the error incured by this approximation is estimated by computing the norm of the energy-momentum tensor. This estimation takes into account the relativistic nature of the matter field and gauges the error in the evolution at the level of the Einstein equations. 

Since GRBoondi is very similar to the GRChombo code, simulations performed with GRBoondi can easily be ported to GRChombo, should full NR simulations be required. This is particularly useful if the backreaction is found to be significant at some point in the simulation. The data from a GRBoondi simulation can be used as initial data for a GRChombo simulation, potentially yielding much better initial data for the full NR simulation.

# Key features of GRBoondi

- **Ease of use**: A central pillar of GRBoondi is its ease of use, relative to other NR software. Many of the basic boiler-plate code and complexities are kept within the source code, allowing the user to have only a basic understanding of the code in order to start researching their problem.

- **Arbitrary choice of background spacetime**: The main parts of GRBoondi are extremely modular, allowing for any arbitrary background to be plugged in, even a numerically computed one. This allows for massive versatility in the source code. Users can swap in and out backgrounds with ease. GRBoondi comes pre-equipped with four background classes ready for use, along with testing suites to verify convergence of each class. These include (starred backgrounds are inherited from GRDzhadzha and are immediately usable in GRBoondi):
  - Minkowski Space
  - Kerr-de Sitter black hole
  - *Boosted Schwarzschild black hole
  - *Kerr black hole

- **Arbitrary modifications to base Proca**: GRBoondi uses specific coding idioms that allow arbitrary modifications to the hard-coded equations of motion. This allows for numerical simulations of any generalized Proca theory. The base model is the standard electromagnetic model, subject to the usual Proca constraint. Any generalized Proca theory can be added on top of this by simply adding the additional pieces from the generalized Proca Lagrangian. Examples of simulating standard Proca and a simple non-linear model [@coates2022intrinsic;@Clough_2022;@_nl_t_rk_2023] are included in the codes repository.

- **Accuracy**:  The metric values and their derivatives are computed exactly at each point, while the matter variables are evolved using 4th-order Runge-Kutta time integration and their spatial derivatives calculated using the same finite difference stencils as GRChombo (6th-order stencils are also available from the GRChombo dependency).

- **Various boundary conditions**: Since GRBoondi inherits from GRChombo, all the boundary conditions in GRChombo are also straight-forwardly applicable to GRBoondi. GRBoondi thus supports:
    - Periodic boundary conditions
    - Sommerfeld out-going wave conditions
    - Static boundary conditions
    - Extrapolation at linear and zeroth order
    - Reflective symmetry (e.g. simulating only the top-half of the computational box due to the background symmetry)

- **Checkpointing**: Since GRBoondi inherits from Chombo, it inherits the checkpointing feature. Simulations can be restarted from a checkpoint file saved during the previous run, allowing for arbitrarily long simulations, even on restricted computing clusters. Moreover, should a simulation fail for any reason, it can be restarted using a previous checkpoint.

- **Diagnostics**: GRBoondi comes equipped with many different diagnostic quantities that can be calculated during the course of a simulation and additional tools for users to calculate their own quantities.
    - *Computation of conserved quantities*. Computation of conserved quantities, such as energy and angular momentum densities and their associated fluxes across a surface can be toggled. GRBoondi can also calculate the energy-momentum tensor trace and its square, for studying the importance of backreaction. Some of the included diagnostic quantities are the square of the Proca field,  conserved energy density, Eulerian energy density, conserved angular momentum, conserved energy density flux, conserved angular momentum flux, trace of the energy-momentum tensor, and square of the energy-momentum tensor. The conserved energy and angular momentum correspond to approximate time and rotational killing vectors. These can be disabled if the spacetime does not possess these. The conserved momentums and fluxes computations are inherited from GRDzhadzha.
    - *Excision of diagnostics in user-specified regions*.  Currently, only spherical excision is implemented, where the user specifies the minimum and maximum radius for the diagnostic quantities. The user can specify which quantities they wish to excise as well, as specified in a parameter file. Generalization to more complicated diagnostic regions is case dependent and left up to the user.
    - *Integration of quantities across the computational grid*. The user can specify which quantities they wish to integrate.
    - *Integration of quantities across spherical surfaces*. The user can specify which quantities they wish to integrate over a spherical surface and the radii of those surfaces as well.

- **Tailored post-processing tools**: GRBoondi comes equipped with custom-built post-processing scripts, leveraging VisIt's python interface. Any of the diagnostic quantities computed during the evolution can be plotted using these tools. It also contains a simple integral plotter, for quickly visualizing the various integrals computed during a simulation, leveraging python's Matplotlib package [@Hunter:2007]

![Some examples of the usage of GRBoondi. The left image is the energy density of a Proca cloud superradiantly excited around a rapidly spinning black hole, superimposed with streamlines of the spatial Proca vector field. The right image is a plot in the xz-plane of the square of the Proca 4-vector in the same background. The spin axis of the black hole is orientated along the z-axis. Superimposed on both images is a slice of the computational grid, clearly showing the refinement hierarchy.](Figures/Combined.png)


# Acknowledgements
LH is supported by funding from the European Research Council (ERC) under the European Unions Horizon 2020 research and innovation programme grant agreement No 801781. LH further acknowledges
support from the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany’s Excellence Strategy EXC 2181/1 - 390900948 (the Heidelberg STRUCTURES Excellence Cluster).

SF is indebted to Katy Clough for her stimulating discussions and extensive input to the code. Katy was an invaluable source of aid for the nuances of NR simulations and some of the finer details of GRChombo. 

The simulations performed as part of this release were carried out on the Baden-Württemberg high-performance computing cluster. The authors acknowledge support by the state of Baden-Württemberg through bwHPC.



# References



