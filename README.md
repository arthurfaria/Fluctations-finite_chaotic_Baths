# Thermalization due to a finite and chaotic bath:

This is a software implementation in `C/C++` that focuses on the thermalization process and the Crooks fluctuation theorem for work. The implementation specifically considers the interaction of a Brownian particle with a bath that is both finite and chaotic.

## Repository content

1. QO_poinc.cpp
 	- Poincarè sections for a quartic oscillator. Hamiltonian dynamics is performed by applying fourth-order symplectic integration, 
	this published by Ref.[1].

2. N_chaotic_baths.cpp
 	- System relaxation in contact to a finite and chaotic heat bath. The latter modeled by a quartic oscillator. Hamiltonian dynamics is performed by applying fourth-order symplectic integration, this published by Ref.[1].

3. Forward_FT.cpp
	- Forward protocol is done in a system coupled to a finite and chaotic heat bath. The latter modeled by a quartic oscillator. Hamiltonian dynamics is performed by applying fourth-order symplectic integration, this published by Ref.[1].
	
4. Reverse_FT.cpp
	- Reverse protocol is done in a system coupled to a finite and chaotic heat bath. The latter modeled by a quartic oscillator.	Hamiltonian dynamics is performed by applying fourth-order symplectic integration, this published by Ref.[1].

For further infos and the physical description of the model used throughout the work, take a look in Ref.[2] (arXiv version attached: '2002.04746.pdf').

References:

[1] https://www.sciencedirect.com/science/article/abs/pii/016727899090019L?via%3Dihub

[2] https://iopscience.iop.org/article/10.1088/1751-8121/ab9a78/meta
