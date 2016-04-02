SOP-GPU
=======

The SOP-GPU package, where SOP-GPU stands for the Self Orginized Polymer Model fully implemented on a Graphics Processing Unit (GPU), is a scientific software package designed to perform Langevin Dynamics simulations of the mechanical unfolding/deformation of large biomolecular systems on the experimental subsecond (millisecond-to-second) timescale. The SOP-GPU package utilizes the Cα-atom based coarse-grained description of proteins combined with Langevin Dynamics in overdamped limit.

The package is fully-implemented on GPU using NVIDIA CUDA technology with focus on high-performance, ease of use, and extensibility [[1][Zhmurov2010],[2][Zhmurov2011]]. SOP-GPU provides out-of-the-box capabilities for numerical simulations of nanoindentation experiments, as well as force-ramp and force-clamp protein pulling. One of the features is optinal support for inclusion of hydrodynamics interactions [[3]][Alekseenko2016].

SOP-GPU have been successfully used to model such system as Fibrin(ogen) molecules [[4]][Zhmurov2011a], CCMV capsid [[5][Kononova2013a],[6][Kononova2016]], microtubule protofilament [[7][Kononova2013b],[8][Theisen2012],[9][Theisen2013]], human synaptotagmin 1 [[10][Duan2011]], and muscle anchoring complex [[11][Bodmer2015]].

Documentation
-------------

Latest documentation is available at [ReadTheDocs](http://sop-gpu.readthedocs.org/en/latest/) ([PDF](http://readthedocs.org/projects/sop-gpu/downloads/pdf/latest/)).


Licensing
---------

This software is distributed under GPLv3 or later (see `COPYING`).

If used for scientific publications, please cite [[1]][Zhmurov2010] and [[2]][Zhmurov2011]. If hydrodyamics functionality is used, please also cite [[3]][Alekseenko2016].


Citations
---------

 1. [Zhmurov, A., Dima, R. I., Kholodov, Y., & Barsegov, V. SOP-GPU: Accelerating biomolecular simulations in the centisecond timescale using graphics processors. Proteins 78, 2984–99 (2010)][Zhmurov2010]
 2. [Zhmurov, A., Rybnikov, K., Kholodov, Y., & Barsegov, V. Generation of random numbers on graphics processors: forced indentation in silico of the bacteriophage HK97. J. Phys. Chem. B 115, 5278–88 (2011)][Zhmurov2011]
 3. [Alekseenko, A., Kononova, O., Kholodov, Y., Marx, K.A., & Barsegov, V. SOP-GPU: influence of solvent-induced hydrodynamic interactions on dynamic structural transitions in protein assemblies. J. Comput. Chem., in press (2016)][Alekseenko2016]
 4. [Zhmurov, A., Brown, A.E.X., Litvinov, R.I., Dima, R.I., Weisel, J.W., & Barsegov, V. Mechanism of fibrin(ogen) forced unfolding. Structure 19, 1615–24 (2011)][Zhmurov2011a]
 5. [Kononova, O., Snijder, J., Brasch, M., Cornelissen, J., Dima, R.I., Marx, K.A., Wuite, G.J.L., Roos, W.H., & Barsegov, V. Structural transitions and energy landscape for Cowpea Chlorotic Mottle Virus capsid mechanics from nanomanipulation in vitro and in silico. Biophys. J. 105, 1893–903 (2013)][Kononova2013a]
 6. [Kononova, O., Snijder, J., Kholodov, Y., Marx, K.A., Wuite, G.J.L., Roos, W.H., & Barsegov, V. Fluctuating nonlinear spring model of mechanical deformation of biological particles. PLOS Comput. Biol. 12, e1004729 (2016)][Kononova2016]
 7. [Kononova, O., Kholodov, Y, Theisen, K.E., Marx, K.A., Dima, R.I., Ataullakhanov, F.I., Grishchuk, E.L., & Barsegov, V. Tubulin bond energies and microtubule biomechanics determined from nanoindentation in silico. J. Am. Chem. Soc. 136(49), 17036–45 (2014)][Kononova2013b]
 8. [Theisen, K.E., Zhmurov, A., Newberry, M.E., Barsegov, V., & Dima, R.I. Multiscale modeling of the nanomechanics of microtubule protofilaments. J. Phys. Chem. B 116(29), 8545–8555 (2012).][Theisen2012]
 9. [Theisen, K.E., Desai, N.J., Volski, A.M., & Dima, R.I. Mechanics of severing for large microtubule complexes revealed by coarse-grained simulations. J. Chem. Phys. 139(12), 121926 (2013).][Theisen2013]
 10. [Duan, L., Zhmurov, A., Barsegov, V., & Dima, R.I. Exploring the mechanical stability of the C2 domains in human synaptotagmin 1. J. Phys. Chem. B 115(33), 10133–10146 (2011).][Duan2011]
 11. [Bodmer, N.K., Theisen, K.E., & Dima, R.I. Molecular investigations into the mechanics of a muscle anchoring complex. Biophys. J 108(9), 2322–2332 (2015).][Bodmer2015] 

[Zhmurov2010]: http://dx.doi.org/10.1002/prot.22824
[Zhmurov2011]: http://dx.doi.org/10.1021/jp109079t
[Alekseenko2016]: http://dx.doi.org/10.1002/jcc.24368
[Zhmurov2011a]: http://dx.doi.org/10.1016/j.str.2011.08.013
[Kononova2013a]: http://dx.doi.org/10.1016/j.bpj.2013.08.032
[Kononova2016]: http://dx.doi.org/10.1371/journal.pcbi.1004729
[Kononova2013b]: http://dx.doi.org/10.1021/ja506385p
[Theisen2012]:  http://dx.doi.org/10.1021/jp212608f
[Theisen2013]: http://dx.doi.org/10.1063/1.4819817
[Duan2011]: http://dx.doi.org/10.1021/jp2025945
[Bodmer2015]: http://dx.doi.org/10.1016/j.bpj.2015.03.036

