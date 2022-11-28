%Example of how to run the ABM, using long timesteps (12 hours)
%If you want to use differently sized timesteps, use the ABM_run_variable
%file. Default parameter settings are provided in the paper and can be
%altered in the get_system_parameter file. 

clear all;close all;clc
addpath('./functions/') %Include all functions in path.
addpath('./subroutines_ABM/') %Add the subroutines to the path

simname = "RUN"; %define name of video for saving

[mySystem,cnst] = getSystemParams_long([125,125]);   %get system parameters

%Adjust initial cell numbers (optional)              
mySystem.params.M1cellNo = 0; %starting number of M1 macrophages
mySystem.params.M2cellNo = 0; %Starting number of M2 macrophages

%Additional modeling settings
mySystem.params.Fradius = 15; %Radius of the circle/ellips of the acinus 
mySystem.params.Shape = 2;    %2 = it will be an ellips, 1 = it will be a circle
mySystem.params.Fellips = 15; %Radius of the second parameter defining ellips geometry (y-direction)
mySystem.params.Number = 1;   %Number of acini that you want to plot on the grid 

%run model 
[mySystemnew,TUcells,~,finalImage,TSC,Tm,Tm2,MutCells,BM] = growTumor(mySystem,cnst);

%Plot some results 
%Plot relative number of tumor cells over time
TUcellNo = TUcells(:,1);
TUnoRelative = 1+((TUcellNo-TUcellNo(1))/(TUcellNo(1))); %calculate relative cell number 
time = linspace(0,cnst.nSteps*0.5,cnst.nSteps+1);        %Create a time variable for plotting 
figure()
plot(time,TUnoRelative,'Color','r','LineWidth',3,'DisplayName','PCABM')
xlabel('Time (days)'); 
ylabel('Relative Number of Tumor Cells')
title("Simulated Tumor Cells over Time")
set(gca, 'FontSize', 12)

%Plot percentage of mutated cells (of all luminal cells)
Ratio2=(Tm2./TUcellNo)*100; %Determine percentage of mutated luminal cells
figure()
plot(time,Ratio2,'Color','r','LineWidth',3), hold on
xlabel('Time (days)');
ylabel('Percentage of mutated cells')
title('Mutated cells to normal cells ratio')

%Save output video
writeMyVideo(finalImage,"subroutines_ABM/vid_"+simname,4) %Create video for saving
saveas(gcf,"subroutines_ABM/plot_"+simname+".png") %Save the video using the earlier provided name