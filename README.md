# CollectiveDecitionMakingModel

In this repo, we provide the necessary code used to regenerate the figures in the paper titled 'Collective Decision-Making Over Nonlinear Potential Fields.' Particularly,
- Figures 2(a)-(d) can be generated using:
  - computations/unbaised_undirected.m
  - computations/unbaised_balanced.m
  - computations/plot_graphs.m
- Figures 3(a)-(d) can be generated using:
  - computations/basins_bais_analysis.m
 
The necessary code corresponding to the simulation of both applications in the paper, Kapferer Tailorshop and fish schools, can be found under simulations/tailorshop and simulations/fish_schools, respectively.  
Further, under the simulations/tailorshop, the data-set mat files are included corresponding to the networks describing the interactions in the tailorshop. mat files ending with 'i' correspont to instrumental (directed) network and 's' corresponds to the sociotional (undirected) networks. 1 indicates the instance at the time of a failed strike attempt, and 2 at indicates the instance at the time of a successful strike. Further elaboration on the data-set can be found here: 

https://www.stats.ox.ac.uk/~snijders/siena/kapfdata.htm
