## Analysis and simulation scripts for: "De Novo DNA-based Catch Bonds"
This depository contains the scripts used in the paper "De novo DNA-based Catch Bonds", by Martijn van Galen, Annemarie Bok, Taieesa Peshkovsky, Jasper van der Gucht, Bauke Albada and Joris Sprakel.
Please note that the scripts are annotated with comments, to instruct the user. Python scripts are provided for: 

**I. Data analysis for the rolling adhesion experiments**

**II. oxDNA simulations and data analysis**

**III. rolling adhesion simulation and data analysis**

Here follows an short explanation for each of these:

## I. Data analysis for the rolling adhesion experiments

The experiments_rollingadhesion folder contains a sample dataset of a measurement with one truncated movie of 25 images (experiment_dataset_snippet) and a scripts folder containing the analysis scripts. The truncated dataset serves to illustrate how to structure the experimental data for the analysis scripts; please note that it is not a full experimental dataset of the manuscript, as this would be too large to upload on github. To run the scripts on your own dataset, please specify the path to the dataset you wish to analyse using the "datafolders" parameters at the top of the scripts, and fill in the other data recording parameters you used in your experiments. Please keep in mind that the scripts should be run in the right order: script 2 uses information extracted by script 1, etc.

To analyse the dataset, the following scripts should be run:

	1.particletracking.py: This script reads the raw rolling adhesion movies, performs pre-processing steps, and stores particle trajectories.

	2.rolling_analysis_separateregimes.py: This script takes the trajectories determined by 1.particletracking.py, separates the trajectories into rolling and stopping regimes, and stores the information of the rolling and stopping regimes as csv files.

	3.analyze_velocities_rollingregimes.py: This script takes the rolling regimes identified by 2.rolling_analysis_separateregimes.py, and makes histograms of the rolling velocities. These histograms are stored as csv files. 4.plot_velocities_rollingregimes.py can be used to plot the average velocities over flowrate.

	4.overlay_velocities_rollingregimes.py: This script takes the mean rolling velocities determined by 3.analyze_velocities_rollingregimes.py and plots the overlayed v_roll vs gammadot graphs (fig. 2E of the manuscript).

	4.plot_velocities_rollingregimes.py: This script takes the mean rolling velocities determined by 3.analyze_velocities_rollingregimes.py. This script only plots one v_roll vs gammadot graph at a time, to plots separate graphs (fig. 4 of the manuscript).

	5.plot_stoptimes_stoppingregimes.py: This script takes the stop times stored in by 2.rolling_analysis_separateregimes.py and plots the mean tstop as a function of gamma dot, overlaying the slip bond and catch bond data (fig. 2f of the manuscript)

	6.overlay_froll_SB_CB.py: This script takes the trajectories found by 1.particletracking.py and determines the fraction of rolling particles f_roll. An overlayed plot of f_roll versus shear rate (gamma dot) is then made and the data of the slip bonds and catch bond dataset (latch 7) is overlayed (fig. 2f of the manuscript)

	6.plot_froll_separate.py: This script takes the trajectories found by 1.particletracking.py, and determines the fraction of rolling particles f_roll. This script plots a separate f_roll over shear rate(gamma dot plot) for each data series: (slip bond, catch bond or one of the latch variations) (fig. 4b of the manuscript).

	7.plot_displacementtraces_article.py: This script takes the rolling and stopping regime information stored by 2.rolling_analysis_separateregimes.py and plots a series of displacement traces, (fig. 2D, and supplementary fig. S4 and S5 in the manuscript)

	8.overlay_traces_movies.py: This script loads the particle trajectories analyzed by 1.particletracking.py and overlays them over the raw movie frames to generate the movies in supplementary videos 3 (slip bonds) and 4 (catch bonds latch7). The colorcoded traces are plotted in a 'traces_SV_tr' subfolder of the dataset, such as: '../Data/example_dataset_snippet/dataset/35.0ulmin/analysis/traces_SV_tr/'

	9.plot_velocityhistograms.py: This script takes the information on the rollig regimes extracted by 2.rolling_analysis_separateregimes.py and plots a separate histogram of vroll for every experminent, overlayed with a Gaussian fit (Fig. S6 of the manuscript).


## II. oxDNA simulations and data analysis
### Simulation:
The simulations_oxDNA folder contains a "data" folder with scripts to initialize a simulation series for simulating the catch bond activation "sim_activate_30pN" (for fig S2a) and catch bond deactivation "sim_deactivate_0pN" (fig S2b).
Each of these folders contains a folder labeled "simulationscripts", which can be used to initiate a simulation series. The "generate_measurearchitecture.py" script allows the user to set the simulation parameters. By running generate_measurearchitecture.py with python, a simulation folder containing oxDNA input.txt files is generated.

The simulation series can then be run in parallel on a computational cluster by using the SLURM task manager, which can execute the series by running "sbatch initiate_jobarrayparallel.sh" from the command line. Alternatively, individual simulations can be run separately, by running "oxdna input.txt" from the command line.
The other scripts in the simulationscripts folders contain functions called by "generate_measurearchitecture.py". More detailed information on the simulation parameters is provided in the generate_measurearchitecture.py scripts.

For more information on using the oxDNA package, please visit: https://dna.physics.ox.ac.uk/index.php/Documentation

### Analysis:
The simulations_oxDNA folder contains the "analysisscripts/scripts_analysis/UTILS/" subfolder. In UTILS, all oxDNA simulation analysis scripts can be found.
Please note that these scripts must be run in python2
The scripts include:

	-1.analyze_states_overtime.py: This script takes a simulation trajectory (mstrajectory.dat) generated during the oxDNA simulations, and computes the connection order parameters of the latch sequence, toehold, cross-link and barrier sequence over the course of the simulations. These results are stored as states.csv

	-2.plot_states_overtime_activation.py: This script reads states.csv generated by 1.analyze_states_overtime.py and plots the f_latch and f_toehold order parameters as a function of time (fig. S2A of the manuscript)

	-3.plot_states_overtime_deactivation.py: Similar to 2.plot_states_overtime_activation.py. Generates a plot of f_latch and f_toehold as a function of time for the deactivation simulation (fig. S2B of the manuscript)

	-4.traj2blender.py: This script is for visualization. it converts an oxDNA simulation trajectory (mstrajectory.dat) into a csv file that can be read by 5.blendersnapshot.py and 6.blender_video.py, to form blender files for vizualisation.

	-5.blendersnapshot.py: This script takes the mstrajectory.dat.csv file generated by 4.traj2blender.py, takes a single snapshot of the simulation and generates a blender workspace for making renders, such as in the insets of fig. S2A and S2B of the manuscript. The user may indicate at the top of the script which snapshot should be stored. This script must be run directly via blender in the command line, as follows: "blender -b -P blendersnapshot.py"

	-6.blender_video.py: Similar to 5.blendersnapshot.py, but this script generates a blender animation, storing a range of simulation snapshots to generate supplementary videos 1 and 2 of the manuscript. The user may indicate at the top of the script which snapshots should be stored. This script must be run directly via blender in the command line, as follows: "blender -b -P blendersnapshot.py"

All other scripts in this folder, including the src folder are auxillary scripts belonging to the oxDNA package that are required for the analysis scripts to function.



## III. rolling adhesion simulation and data analysis

The simulations_rollingadhesion folder contains the following subfolders:

**catchbond_heatmap_simulations:** Contains all scripts required to perform the simulation series to generate the catch bond activation heat map (fig. 3B)

**catchbond_rollsimulations** Contains all scripts required to perform the simulatio nseries to generate the catch bond v_roll and t_stop plots (fig. 3C,D,E)

**slipbond_rollsimulations** Contains all scripts required to perform the simulation series to generate the slip bond v_roll and t_stop plots (fig. 3C,D,E)

Each of these simulation series contains the following subfolders:

Simulation_series_name/
	simulationscripts/
	analysis/

All scripts required to run the simulations are stored in simulationscripts/
The analysisscripts/ folder contains all the scripts required to analyse the raw data and processed results after running the analysisscripts will be stored in the analysis/ folder.

Here, we will go over all the simulation scripts in the simulationscripts/ folder required to a. generate the input parameters for a simulation series, b. run the simulations, and c. analyze the data

### a. Generate input parameters:
-generate_measurearchitecture.py: This file allows the user to specify all the simulation parameters, and generates a number of input files that can be read by the run_simulation.py script. In the top of the file, a list of variable parameters is specified. Below that, the remaining constant parameters are specified. More detailed information on the parameters is provided as comments in the generate_measurearchitecture.py file. When run using "python generate_measurearchitecture.py", this script creates the ../simulation/ folder and generates an input.txt file containing the parameters for each individual simulation.
Furthermore, a simulationslist.txt file is generated, which is a list of all the simulations that need to be performed. 

-generate_input.py: This script is called by generate_measurementarchitecture.py

-generate_simulationslist.py: This script is called by generate_measurearchitecture.py

### b. Running the simulations:
The main simulation script is called run_simulation.py. This script is used to start a simulation. An individual simulation can be run with this script, using "python run_simulation.py ../simulation/path_to_simulationfolder/ input.txt"

-rollingbead_f.py:
This is the main rolling adhesion simulation function, which is called by run_simulation.py. An explanation of all the steps is given in the script. rollingbead_f.py calls several other scripts that perform part of the simulation:

	-InitializeArray.py: 		This function creates initial positions for the linkers on the channel surface. These positions are chosen such that the bond density rho is met. Next, the function randomly allows some linkers to form
 	-ShearForceAndTorque.py: 	This function computes the Force and Torque on the particle caused by the fluid shear.
  	-Equilibrate.py:		This function equilibrates the free energy of the system by rotating and translating the bead using a grandient descent approach
   	-RotateBead.py: 		This function is called by Equilibrate.py. It is used to recalculate the positions p1 of the linkers on the particle as the particle rotates.
    -TotalLinkerForceAndTorque.py:	This function is called by both Equilibrate.py and rollingbead_f.py. It computes the total extensional force F and torque T all bound linkers apply to the particle, the relative extension x, and finaly the x component of the force F. 
    -UpdateArray.py:		This function is called by rollingbead_f.py at the start of each simulation loop. It creates and removes positions of linkers on the channel surface in case the bead has moved during the previous step. Parts of UpdateArray.py are similar to InitializeArray.py
    -RateConstants.py:		This function is called by rollingbead_f.py before each kinetic monte carlo step, to compute the rate constants for formation and dissociation of every linkers.

These scripts vary slightly between the catch bond and slip bond simulation folders, as catch bond simulations allow bonds to interconvert between states, whereas bonds in the slip bond simulations can only exist in one state.

The batch of simulations generated by generate_measurearchitecture.py can be run using the SLURM workload manager "https://slurm.schedmd.com/documentation.html", by executing the initiate_jobarrayparallel.sh file.
This is done using the command line: "sbatch initiate_jobarrayparallel.sh simulationslist.txt"
The initiate_jobarrayparallel.sh script reads in the list of simulations specified in simulationslist.txt, and then runs all simulations. Please make sure a sufficient number of tasks (--array command) is specified at the top of the initiate_jobarrayparallel.sh script. Alternatively, individual simulations can be run separately by running python from the command line, with the command: "python run_simulation.py ../simulation/dil_X/stressX/repX/ input.txt".

### c. Data analysis
The following data analysis scripts are found in the analysisscripts subfolders:

#### catchbond_heatmap_simulations/analysisscripts/:
	-plot_heatmap_activationdistribution.py: This script plots the average heatmap of the distribution of weak and strong bonds over the course of a rolling adhesion catch bond simulation (fig. 3B of the manuscript)

#### catchbond_rollsimulations/analysisscripts/ and slipbond_rollsimulations/analysisscripts/:
Please run these in order.

	-1.stop_versus_roll_interp.py: 		This script takes the particle location trajectory over the simulation series and separates it into stopping and rolling regimes. As the rolling particle simulations are event-driven, they do not provide location information that is linearly spaced in time. The script therefore performs a linear interpolation of the location information prior to regime separation.
						This script can also be used to plot displacement traces (figs. S8 and S9.)

	-2.analyze_rollingvelocities.py:	This scripts takes in rolling regimes analyzed by 1.stop_versus_roll_interp.py and computes the weighted gaussian mean and standard deviation of the v_roll in the rolling regimes.

	-3.overlay_vroll_vs_shearrate.py:	This scripts takes in the average rolling velocities analysed by 2.analyze_rollingregimes.py for both the slip bond simulation and the catch bond simulation, and plots v_roll vs shearrate for both of these (fig. 3D of the manuscript). Please make sure the simulations of both the slip bond dataset and the catch bond dataset have been fully analysed by 1.stop_versus_roll_interp.py, followed by 2.analyze_rollingregimes before running this script.

	-4.analyze_tstop.py:			This scripts takes in rolling regimes analyzed by 1.stop_versus_roll_interp.py and computes the average tstop and the standard deviation of tstop.
	
	-5.overlay_tstop_vs_shearrate.py:	This scripts takes in the average tstops analysed by 5.analyze_tstop.py for both the slip bond simulation and the catch bond simulation, and plots t_roll vs shearrate for both of these (figure fig. 3E of the manuscript).




