SOP-GPU
=======

SOP-GPU software package is used for coarse-grained simulations of proteins and other biological macromolecular assemblies.

It is fully-implemented on GPU using NVIDIA CUDA technology, with focus on high-performance, ease of use, and extensibility [[1][Zhmurov2010],[2][Zhmurov2011]]. SOP-GPU provides out-of-the-box caabilities for numerical simulations of nanoindentation experiments, as well as force-ramp and force-clamp protein pulling. One of the features is optinal support for inclusion of hydrodynamics interactions [[3]][Alekseenko2016].

SOP-GPU have been successfully used to model such system as Fibrin(ogen) molecules [[4]][Zhmurov2011a], CCMV capsid [[5][Kononova2013],[6][Kononova2016], and microtubule protofilament [[7]][Theisen2013].

Documentation
-------------

Latest documentation is available at [ReadTheDocs](http://sop-gpu.readthedocs.org/en/latest/) ([PDF](http://readthedocs.org/projects/sop-gpu/downloads/pdf/latest/)).


Licensing
---------

This software is distributed under GPLv3 or later (see `COPYING`).

If used for scientific publications, please cite [[1]][Zhmurov2010] and [[2]][Zhmurov2011]. If hydrodyamics functionality is used, please also cite [[3][Alekseenko2016]]


Citations
---------

 1. [Zhmurov, A., Dima, R. I., Kholodov, Y. & Barsegov, V. SOP-GPU: accelerating biomolecular simulations in the centisecond timescale using graphics processors. Proteins 78, 2984–99 (2010)][Zhmurov2010]
 2. [Zhmurov, A., Rybnikov, K., Kholodov, Y. & Barsegov, V. Generation of random numbers on graphics processors: forced indentation in silico of the bacteriophage HK97. J. Phys. Chem. B 115, 5278–88 (2011)][Zhmurov2011]
 3. [Alekseenko, A., Kononova, O., Kholodov, Y., Marx, K.A. & Barsegov, V. SOP-GPU: influence of solvent-induced hydrodynamic interactions on dynamic structural transitions in protein assemblies. J. Comput. Chem., in press (2016)][Alekseenko2016]
 4. [Zhmurov, A. et al. Mechanism of fibrin(ogen) forced unfolding. Structure 19, 1615–24 (2011)][Zhmurov2011a]
 5. [Kononova, O. et al. Structural transitions and energy landscape for Cowpea Chlorotic Mottle Virus capsid mechanics from nanomanipulation in vitro and in silico. Biophys. J. 105, 1893–1903 (2013)][Kononova2013]
 6. [Kononova, O. et al. Fluctuating Nonlinear Spring Model of Mechanical Deformation of Biological Particles. PLOS Comput. Biol. 12, e1004729 (2016)][Kononova2016]
 7. [Theisen, K. E., Desai, N. J., Volski, A. M. & Dima, R. I. Mechanics of severing for large microtubule complexes revealed by coarse-grained simulations. J. Chem. Phys. 139, 121926 (2013)][Theisen2013]

[Zhmurov2010]: http://dx.doi.org/10.1002/prot.22824
[Zhmurov2011]: http://dx.doi.org/10.1021/jp109079t
[Alekseenko2016]: http://dx.doi.org/10.1002/jcc.24368
[Zhmurov2011a]: http://dx.doi.org/10.1016/j.str.2011.08.013
[Kononova2013]: http://dx.doi.org/10.1016/j.bpj.2013.08.032
[Kononova2016]: http://dx.doi.org/10.1371/journal.pcbi.1004729
[Theisen2013]: http://dx.doi.org/10.1063/1.4819817
