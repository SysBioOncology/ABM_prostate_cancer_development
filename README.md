# ABM_prostate_cancer_development
Code used in the paper "Multicellular agent-based model of prostate tumor onset and development" by Passier et al.

### General description 
In this study, we created a 2D Agent-Based Model (ABM) to simulate prostate cancer onset and development in a prostatic acinus. The model includes the following agents: basement membrane, basal cells, fibroblasts, luminal cells (capable of gaining mutations), stem cells, M1 and M2 macrophages, Cancer Associated Fibroblasts and stroma. All in silico experiments were conducted using a 125x125 grid and timesteps of 12 hours. Model parameters were fitted to experimental in vitro data of cell cultures of the tumor cells (luminal cells), fibroblasts, M1 and M2 macrophages. Parameters for the other cell types were based on literature (for a detailed description of the derivation of all model parameters, please see the supplementary tables of the paper). 

To validate the model simulations, in silico growth patterns were compared to growth patterns often observed in patients (in collaboration with a pathologist and oncologist). Additionally, we created eight different patient phenotypes in a dataset of prostate cancer patients from the TCGA. We defined three patient markers ( mutational burden, macrophage tumor stimulation and CAF tumor stimulation) and classified the patients as 'high' or 'low' for all three patient markers (based on measured quantities: e.g. macrophage tumor stimulation was defined by: amount of M2 macrophages, level of CXCL2 and level of STAT3; chemokines involved in the transition towards a pro-tumor phenotype). This yields eight classes, ranging from 'high' for all three markers to 'low' for all three markers. The corresponding model parameters were adapted as well, so that these same eight phenotypes could be simulated with the model. The progression free survival times for these eight groups of prostate cancer patients were compared to the predicted tumor loads for the simulated patients after invasion (so when the simulated patient would be classified as having cancer). We found a significant correlation between predicted tumor load and actual disease progression in patients.  

### Uses of the model 
The code can be used to run simulations of the onset and development of prostate cancer in a prostatic acinus. Model parameters can be adjusted to mimick specific classes of patients and a multitude of outputs can be plotted and studied. In order to run the general model, the function *Example_run_ABM* needs to be in the main folder. Additionally, the user needs to place the 'subroutines_ABM' folder in this main folder. The *Example_run_ABM* script can then be run: it will show the video of the simulation and plot the amount of tumor cells and mutations over time upon finishing the simulation. Should a different time step be desired, the *ABM_run_variable* script can be used. The desired timestep can be provided in this script. Lastly, to run simulations for all described eight phenotypes, the file *Eight_phenotypes_run* should be used. This file runs the script and plots average behavior (including standard deviation) for all eight phenotypes

### Software 
The code was developed and simulations were conducted using Matlab 2021a. For full usage of the code, the following packages should be installed: the Image Processing Toolbox, the Statistics and Machine Learning Toolbox, the optimization toolbox and the Global Optimization Toolbox (for fitting). 

### Files 

### folder 'subroutines_ABM'
Contains all files to run the prostate cancer agent based model (PCABM). 

- **getSystemParams_long**: Retrieves all general parameters needed for model simulations, for simulations with timesteps of 12 hours. If you want to change a parameter, you can do it here. 
- **getSystemParams_variable**: Retrieves all general parameters needed for the model simulations, adjusting them to any timestep (lower than 15 hours). You should provide this timestep in ABM_run_variable

- **growTumor**: The main function of the ABM, that bundles all other functions and runs all parts. It defines the starting geometry and then loops through all the timesteps, allowing for all agents to take action each iteration. 
- **growTumor_variable**: The main function if you use variable timesteps.  

- **B_die**: Basement membrane cannot die or proliferate on its own, but a round for this celltype was created to allow the user to adapt these properties for future use (or if you want to prescribe any other actions to the basement membrane). 

-  **C_go_grow_die**:Script for basal cells death and proliferation.

- **F_go_grow_die**: In principle, fibroblasts are modelled as quiescent cells. However, they could proliferate, migrate and die if the user wishes to implement this. 
Furthermore, this is the round in which fibroblasts can transition towards CAFs. 

- **CAF_go_grow_die**: Function that models proliferation, migration and death of CAFs. Additionally, CAFs can elicit EMT in tumor cells during this round (enabling mutated cell migration) and break down ECM. Lastly, they can affect macrophage differentiation towards M2 macrophages when in range and promote tumor cell proliferation (with the function CF_promote). 

- **CF_promote**: Function that models tumor cell proliferation promotion by CAFs. 

- **getAdjacent**: Get cells in the Moore neighborhood of the current cell (to determine whether cells can migrate or proliferate towards adjacent spaces) and create masks; random number vectors (both for all agents and for agents that have a free spot in their neighborhood)

- **getAdjacent2**: Get masks for cells in a two times Moore neighborhood. (neighbors for all cells and random number vectors). 

- **mCellRound**: Main program for the M1 macrophages, this calls on all the other macrophage functions (so that all macrophage actions are performed during this round: migration, death and killing of the tumor cells).

- **M_go2**: Function that allows for migration by macrophages of more than one gridspace (up to 2 gridspaces away) at a time (to enable invasion of the acinus). 

- **M_kill**: Function for modelling mutated cell killing by macrophages (both M1 and M2). 

- **M_go_die**: Function that models the final action for the macrophages at the end of their round (they can migrate more than once, which is done using the M_go function, die or idle). 

- **mCellRound2**: Main program for the M2 macrophages, also calls all the other macrophage functions, for the M2 macrophages (allowing M2 macrophages to die, migrate, kill, promote TU cell proliferation and stimulate EMT of tumor cells). 

- **M_promote**: Function for modelling mutated cell proliferation promotion by M2 macrophages. 

- **mutation_round**: Function that models luminal cells gaining mutations and the proliferative advantage that comes with the mutations (it also takes into account the increased chance of gaining a mutation for cells that have already mutated. 
- **mutation_round_variable**: Similar function, but then for the variable time step. 

- **shuffleCells**: Randomize the order vector of cells (so that you start with different cells each time, do not have the same order each round)

- **TU_go_grow_die**: Function for modeling tumor cell migration, growth and death.
Additionally the function to break down basement membrane is called during this function. 

- **TU_kill**: Function that models basement membrane break down by mutated luminal cells. 

- **TU_interaction_round**: Function that models break down of the ECM by mutated luminal cells and affected macrophage differentiation. 

- **updateSystem**: Function to update sytem in each time step (iteration) after all cells performed their round of actions.

- **visualizeSystem**: Visualization of the simulation. 
- **visualizeSystem_variable**: Same function, but then using a variable time step

- **writeMyVideo**: Function to write visualized simulation to video file

- **Ploterr**: Function to visualize average behavior of a number of simulations and shaded error of the mean, standard deviation or distribution. 
Written by Brendan Hasz, retrieved from: https://github.com/brendanhasz/matlab-uncertainty-viz/blob/master/ploterr.m 

The second folder, called 'functions' contains the files used for fitting the model parameters to the data. It shows the optimization programs and script to run the Particle Swarm Optimization algorithm. Beware that a data folder and the data file need to be added before this script can be run. 

### folder 'functions'
Contains functions for optimizing parameters (for the functions needed to run the fitting algorithm used for optimizing model parameters; fitting them to the data). 
- **ABM_PSO_function**: Function that runs particle swarm optimization (PSO). PSO fits parameters by comparing best Mean Squared Error (MSE) between Relative Tumor Cell Numbers from data to model output for each fit.
- **MSE**: Calculating MSE for the fitting procedure PSO - determines the difference between the model output for fitted parameter values and actual value in the data and returns this value.
