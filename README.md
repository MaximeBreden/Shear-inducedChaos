This repository contains the Matlab codes used for the paper "Computer-assisted proof of shear-induced chaos in stochastically perturbed Hopf systems", 
by M. Breden and M. Engel. 

All the computer-assisted parts of the proof discussed in the paper can be reproduced by running script_proof.m. The proofs make crucial usage of
interval arithmetic via the Intlab toolbox (http://www.ti3.tu-harburg.de/intlab/). The curious reader can also experiment with the non-validated 
computations perfomed in script_numerics.m, which provides curves describing the influence of the noise level sigma on the Lyapunov exponent, as 
discussed in section 4.4 of the paper.
