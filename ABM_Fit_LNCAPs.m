%Example of how to run parameter optimiziations (Fit them to data using Particle Swarm Optimization)
%Data is structured in an excel file, using three columns: one containing
%the time of measurement (t14 in this case), one containing the amount of
%tumor cells (m14 in this case) and one containing the standard deviation
%of each measurement (s14 for this case). 

%clear all;close all;clc
addpath('./functions/') %Include all functions in path.
addpath('./data/')      %Add the folder in which the data is stored to the path
addpath('./subroutines_ABM/') %Add the folder in which the subroutines are saved to the path

%%
%load data and set some variables
data=readtable('LNCaP_R1881.xlsx','Sheet','LNCaP'); %Specifying name of the file, if the excel file contains multiple sheets, you can indicate the wanted sheet by using this code. 
gridsize = [125,125];               % Standard used grid size for all model simulations
M1cellNo = 0;                       % For this optimization, only tumor cells (LNCAP) were included. If there are also macrophages or fibroblasts, their initial values should be increased. 
M2cellNo = 0;
FcellNo = 0;
varNames = ["TUpprol","TUpmax"];    %select which parameters to optimize
lb = [0,0];                         %set lower fitting bounds for these parameters
ub = [0.5,10];                      %set upper fitting bounds for these parameters
version = "TU";                     %optimizations for tumor cell number ("TU") or fibroblast cell number ("F") - We also conducted a fit containing just fibroblasts. Here you can indicate which scenario you are using. 
optno=50;                           %number of optimization runs
expname='R1881';                    %parameter condition (Here you can specify the data you are using to fit the parameters). 

%% 
t14 = data.t14(~isnan(data.t14));  %Read in timesteps of each measurement
m14 = data.m14(~isnan(data.m14));  %Read in measurements (relative number of tumor cells at each timestep)
s14 = data.s14(~isnan(data.s14));  %Read in standard deviation for each time step
nSteps14 = 34;                     %Number of time steps you want to run the simulation for (adjust this on the basis of your measurements
TUcellNo14 = 1500;                 %Initial number of tumor cells 
name14 = "LNCaP_R1881_Inc14";      %Specify the name of the fit

[mySystem14,cnst14,valVars14,fvals14]=ABM_PSO_function(t14,m14,gridsize,nSteps14,varNames,lb,ub,TUcellNo14,M1cellNo,M2cellNo,FcellNo,version,name14,optno);
%The ABM_PSO_function (in the folder 'functions') details the settings for
%the Particle Swarm Optimization (PSO). 

varsmed=median(valVars14); %Determining the median fitted values 
%fvals indicates the error per fitted value 

