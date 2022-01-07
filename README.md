# Recoverability_of_ARG_topologies
Matlab code to accompany the paper of the same name.

The code for each problem is most often divided into 3 sections: 
-(1) a return function will return the value of a probability matrix entry if feasible indices are called, or 0 otherwise. This is a useful helpful function which keeps the main code tidy,
-(2) a solve function iterates through a given series of recursions and returns a solution matrix,
-(3) a plotter function extracts data from the solution matrix and was used to create all the graphs used in the paper.

Solving the topology of the full ARG including galled recombinations requires some additional code used for optimising an otherwise unfeasibly large matrix. A _reference matrix_ is used to associate the tuple of states for the left (or right) hand lineages of the graph to an integer. The exploits the constraint that the total number of left (or right) hand lineages must sum to the number of open recombination loops. This folder therefore additionally includes:
-(4) creation of the reference matrix,
-(5) return function which return the associated integer to a tuple of states.
