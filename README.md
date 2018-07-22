# Adiabatic Quantum Annealling Simulation

This project is the base code I used to complete my masters' dissertation on the effects of the scheduling on adiabatic annealing in 2016-17.

At the heart of this code is the `Anneal` class which is used to define and simulate a quntum annealing proces in the absence of environmental noise. Into this object you can feed the annealing schedule you would like to explore as well as the connectivitiy between qubits which define the computation. The anneal object is designed to give the user an end to end object, with the ability to build the Hamiltonian for the problem given the interaction coefficients, to run the anneal with any defined schedule input, and to finally plot the outcome of the anneal.

