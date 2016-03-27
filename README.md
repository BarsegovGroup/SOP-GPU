SOP-GPU
=======

SOP-GPU software package is used for coarse-grained simulations of proteins and other biological macromolecular assemblies.

It is fully-implemented on GPU using NVIDIA CUDA technology, with focus on high-performance, ease of use, and extensibility [[1][Zhmurov2010],[2][Zhmurov2011]]. SOP-GPU provides out-of-the-box caabilities for numerical simulations of nanoindentation experiments, as well as force-ramp and force-clamp protein pulling. One of the features is optinal support for inclusion of hydrodynamics interactions [[3]][Alekseenko2016].

SOP-GPU have been successfully used to model such system as Fibrin(ogen) molecules [[4]][Zhmurov2011a], CCMV capsid [[5][Kononova2013a],[6][Kononova2016]], and microtubule protofilament [[7]][Kononova2013b].

Documentation
-------------

Latest documentation is available at [ReadTheDocs](http://sop-gpu.readthedocs.org/en/latest/) ([PDF](http://readthedocs.org/projects/sop-gpu/downloads/pdf/latest/)).


Licensing
---------

This software is distributed under GPLv3 or later (see `COPYING`).

If used for scientific publications, please cite [[1]][Zhmurov2010] and [[2]][Zhmurov2011]. If hydrodyamics functionality is used, please also cite [[3]][Alekseenko2016].


Citations
---------

 1. [Zhmurov, A., Dima, R. I., Kholodov, Y. & Barsegov, V. SOP-GPU: Accelerating biomolecular simulations in the centisecond timescale using graphics processors. Proteins 78, 2984–99 (2010)][Zhmurov2010]
 2. [Zhmurov, A., Rybnikov, K., Kholodov, Y. & Barsegov, V. Generation of random numbers on graphics processors: forced indentation in silico of the bacteriophage HK97. J. Phys. Chem. B 115, 5278–88 (2011)][Zhmurov2011]
 3. [Alekseenko, A., Kononova, O., Kholodov, Y., Marx, K.A. & Barsegov, V. SOP-GPU: influence of solvent-induced hydrodynamic interactions on dynamic structural transitions in protein assemblies. J. Comput. Chem., in press (2016)][Alekseenko2016]
 4. [Zhmurov, A., Brown, A.E.X., Litvinov, R.I., Dima, R.I., Weisel, J.W., & Barsegov, V. Mechanism of fibrin(ogen) forced unfolding. Structure 19, 1615–24 (2011)][Zhmurov2011a]
 5. [Kononova, O., Snijder, J., Brasch, M., Cornelissen, J., Dima, R.I., Marx, K.A., Wuite, G.J.L., Roos, W.H., and Barsegov, V. Structural transitions and energy landscape for Cowpea Chlorotic Mottle Virus capsid mechanics from nanomanipulation in vitro and in silico. Biophys. J. 105, 1893–903 (2013)][Kononova2013a]
 6. [Kononova, O., Snijder, J., Kholodov, Y., Marx, K.A., Wuite, G.J.L., Roos, W.H. & Barsegov, V. Fluctuating nonlinear spring model of mechanical deformation of biological particles. PLOS Comput. Biol. 12, e1004729 (2016)][Kononova2016]
 7. [Kononova, O., Kholodov, Y, Theisen, K.E., Marx, K.A., Dima, R.I., Ataullakhanov, F.I., Grishchuk, E.L., & Barsegov, V. Tubulin bond energies and microtubule biomechanics determined from nanoindentation in silico. J. Am. Chem. Soc. 136(49), 17036–45 (2014)][Kononova2013b]

[Zhmurov2010]: http://dx.doi.org/10.1002/prot.22824
[Zhmurov2011]: http://dx.doi.org/10.1021/jp109079t
[Alekseenko2016]: http://dx.doi.org/10.1002/jcc.24368
[Zhmurov2011a]: http://dx.doi.org/10.1016/j.str.2011.08.013
[Kononova2013a]: http://dx.doi.org/10.1016/j.bpj.2013.08.032
[Kononova2016]: http://dx.doi.org/10.1371/journal.pcbi.1004729
[Kononova2013b]: http://dx.doi.org/10.1021/ja506385p
