function [mySystem,cnst,valVars,fvals]=ABM_PSO_function(t,m,gridsize,nSteps,varNames,lb,ub,TUcellNo,M1cellNo,M2cellNo,FcellNo,version,name,optno)
%This function is called by ABM_Fit_LNCAPs and is used for fitting
%parameters. Here the settings for Particle Swarm Optimization are defined
%and the program loops through the amount of optimizations you want to
%conduct. You should provide data to fit the parameters to and initial
%settings for the fit and parameters to be fitted, and the function returns
%fitted values and corresponding errors of the fit. 

%Make sure everything is added to the path
addpath('./functions/') %Include all functions in path.
addpath('./data/') %Add data to the path
addpath('./subroutines_ABM/') %Add the subroutines to the path

opt = true;                                    %set if the optimization should run several times                                  

[mySystem,cnst] = getSystemParams_variable(gridsize);   %Optain all system parameters
%It is important to choose a timestep that fits the dataset you have, OR
%adjust the parameters according to the timestep after fitting. 

cnst.nSteps = nSteps;                           %set the number of steps for the simulation 
cnst.video = false;                             %don't show the video on screen (delays the simulations
cnst.simnumber = 1;

%Save provided parameters for this function 
mySystem.params.TUcellNo = TUcellNo;            %set number of initial tumor cells
mySystem.params.M1cellNo = M1cellNo;            %number of initial M1 cells
mySystem.params.M1influxProb = 0;               %no macrophage influx (depending on the simulation: in the example there is only tumor cells, so we do not allow influx of macrophages)
mySystem.params.M2cellNo = M2cellNo;            %number of initial M1 cells
mySystem.params.FcellNo = FcellNo;              %number of initial fibroblasts 
mySystem.params.TUps = 0;                       %no symmetric division as there are no stem cells in vitro (depends on the used data)

costFunction = @(x) MSE(x,mySystem, cnst,  varNames, version, t, m);   %create cost function for PSO
%Provides a function that allows for calculation of the error between the fitted
%simulation and data so that the algorithm can minimize it. 

options = optimoptions('particleswarm','Display','iter');       %set some options for the PSO

%It loops through the number of optnumbers provided 
if opt       
    for i=1:optno
        [vals,fval,exitflag,output] = particleswarm(costFunction, length(varNames), lb, ub,options); %PSO
        fvals(i)=fval; %Error for this optimization
        valVars(i,:)=vals; %Fitted value(s) for this optimization 
        disp("Finished simulation "+name+" "+num2str(i))
        save (name+"_valVars.mat",'valVars') %Store the found fitted values
        save (name+"_fvals.mat",'fvals') %Store the found error values, relating to the fitted values 
    end
else
    [valVars,fval,exitflag,output] = particleswarm(costFunction, length(varNames), lb, ub,options); %PSO
end

