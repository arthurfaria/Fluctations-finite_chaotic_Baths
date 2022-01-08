Brief files description:

1. QO_poinc.cpp
 	- Poincarè sections for a quartic oscillator. Hamiltonian dynamics is performed by applying fourth-order symplectic integration, 
  this published by Ref.[1].

2. N_chaotic_baths.cpp
 	- System relaxation in contact to a finite and chaotic heat bath. The former modeled by a a quartic oscillator.
  Hamiltonian dynamics is performed by applying fourth-order symplectic integration, this published by Ref.[1].

3. Forward_FT.cpp
	- stochastic dynamics for a Brownian particle under a harmonic potential whose center of mass is displaced with constant velocity (see Ref.[2]). 
  Fluctuation Theorem for work and heat are calculated. Algorithm based on Ref.[3]
	
4. Reverse_FT.cpp
	- generalized Fokker-Planck (GenBM) equation (see 'gener_FP.pdf' file for further infos). Both the generalized semiclassical distriubtion (GenBM_rho) 
  and the distribution of a standard Brownian motion (BM_rho) are obatined considering a external harmonic potential. 
  Algorithm based on finite diference approach to compute derivatives.

For further infos and the physical description of the model used throughout the work, take a look in Ref.[2] (arXiv version attached).

References:
[1] https://www.sciencedirect.com/science/article/abs/pii/016727899090019L?via%3Dihub
[2] https://iopscience.iop.org/article/10.1088/1751-8121/ab9a78/meta
