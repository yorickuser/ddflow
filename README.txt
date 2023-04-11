######################################################
###     Program files for numerical analysis of    ###
###     innovation-driven taxonomic turnovers      ###
######################################################
Written by Hiroshi C. Ito (2021)
Email: hiroshibeetle@gmail.com


This directory contains two program files: "ddflow_osm1_amnt.R", "ddflow_psmb0b_amnt.cc".

These programs were used for producing numerical results in a paper titled "The adaptation front equation explains innovation-driven taxonomic turnovers and living fossilization" (see the end of this file for the abstract).

See the description below for execution of these programs.

################################################
##            ddflow_osm1_amnt.R              ##
################################################

<Description>
R program for simulation of DD-flows (innovation-driven taxonomic turnovers) with OSM
(oligomorphic stochastic model, Ito and Dieckmann (2014)).

<Requirement>
R-package "simevol-0.1.4" is required.
simevol requires R-package "deSolve".
simevol can be installed by using devtoos:
library(devtools)
install_github('yorickuser/simevol@0.1.4')

<Execution>
Eco-evolutionary setting can be chosen from the list in "simsets" (e.g., "3_K_sharp_peak" is selected by setid=3).
Execution method 1 (at R prompt): source("ddflow_osm1_amnt.R")
Execution method 2 (at console): Rscript ddflow_osm1_amnt.R setid=3

Simulation can be controlled by updating "command1.R" in an automatically created directory ".simevol/".
For example, from the console in the same directory running "ddflow_osm1_amnt.R",
you can stop the simulation by
  $echo "halt()" > .simevol/command1.R
or exchange the horizontal and vertical axes in the main plotting window by
 $ echo "plot_var(2,1)" > .simevol/command1.R


Values for setid correspond to figures in the paper as follows.
1_default: Fig. 1, 2, 3 
2_default_robust: Fig. S2 
3_K_sharp_peak: Fig. S3
4_K_flat_top: Fig. S4
5_K_bimodal: Fig. S5
6_grady_depend_x: Fig. S10
7_kernel_platykurtic: Fig. S13
8_kernel_leptokurtic: Fig. S12
9_kernel_asymetric: Fig. S14
10_mutation_y_rare: Fig. 4, S18 
11_multiple_geographic_regions: Fig. 5, 6, S19, S20
12_grady_avg: Fig. Fig. S8, S9
13_grady_saturation: Fig. S6, S7

<Operation environment>
OS: Utunbu 18.04
R: R-3.4.4 
R-packages: simevol-0.1.2, deSolve-1.28
(confirmed to work also in Ubuntu 20.04 with R-4.2.0, simevol-0.1.2, and deSolve-1.35)


################################################
##           ddflow_psmb0b_amnt.cc            ##
################################################

<Description>
C++ program for simulation of DD-flows with PSM
(Polymorphic stochastic model, Dieckmann and Law (1996))
for sexual populations (see Ito and Dieckmann (2007)).
The calculated results were used for Figs. S1 and S15.

<Requirement>
X11 library is required for compilation.

<Execution>
Compilation:
g++ ddflow_psm0b_amnt.cc -lX11 -L/usr/X11R6/lib  -lm -lstdc++

Execution (species-level output): ./a.out -u > data.dat
Execution (phenotype-level output): ./a.out -p > data.dat



<Operation environment>
OS: Utunbu 18.04
g++: g++-7.5.0
X11: libx11-dev-2:1.6.4-3ubuntu0.4


################################################
##            Paper information               ##
################################################

Title: The adaptation front equation explains innovation-driven taxonomic turnovers and living fossilization

Author: Hiroshi C. Ito and Akira Sasaki

Abstract:
Evolutionary taxonomic turnovers are often associated with innovations beneficial in various ecological niches. Such innovations can repeatedly occur in species occupying optimum niches for a focal species group, resulting in their repeated diversifications and species flows from optimum to suboptimum niches, at the expense of less innovated ones. By combining species packing theory and adaptive dynamics theory, we develop an equation that allows analytical prediction for such innovation-driven species flows over a niche space of arbitrary dimension under a unimodal carrying capacity distribution. The developed equation and simulated evolution show that central niches (with the highest carrying capacities) tend to attain the fastest innovation speeds to become biodiversity sources. Species that diverge from the central niches outcompete the indigenous species in peripheral niches. The outcompeted species go extinct or evolve directionally toward far more peripheral niches. Due to this globally acting process over niches, species occupying the most peripheral niches are the least innovated and have deep divergence times from their closest relatives, and thus they correspond to living fossils. The extension of this analysis for multiple geographic regions shows that living fossils are also expected in geographically peripheral regions for the focal species group.
