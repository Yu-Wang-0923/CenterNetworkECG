# Multiplex depth for network-valued data


###we provide the main codes for reproducing the simulation results in the paper.

#main_function.R  file provides the main functions for the paper, which should be run first.



#######Simulation codes for unweighted network:

#Gen_unweighted_network.R file provides the functions for generating unweighted networks from SBM, ER,RGG,SW. (should be run before Sim_Unweighted_EA.R and  simulation_averdage rank_unweighted.R)

#Sim_Unweighted_EA.R provides all codes to reproduce the results for finding centers in Figure 2.

#simulation_averdage rank_unweighted.R file  provides the codes to reproduce and plot the results in Figure 3 and Table 3.



#########Simulation codes for weighted network:

#Gen_weighted_network.R file provides the functions for generating weighted networks from models in paper. (should be run before Sim_weighted_EA.R and  average_rank_weighted_network.R)

#Sim_weighted_EA.R file provides all codes to reproduce the results for finding centers in Figure 4.

#average_rank_weighted_network.R file provides the codes to reproduce and plot the results in Table 5.



#########Simulation codes for running time:

#Running_Time.R file provides codes to reproduce the results for Figure 5.



#########Code for plotting results:

#PlotResults.R file provides codes to plot the results for Figures 2,4 and 5.