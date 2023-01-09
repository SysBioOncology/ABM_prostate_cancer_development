%Run the ABM using a timestep of your choosing (be cautious with time steps
%larger than 12 hours, then probabilities will start to exceed 1 and
%results will be distorted. 

%clear all;close all;clc
addpath('./functions/') %Include all functions in path.
addpath('./subroutines_ABM/') %Add the subroutines of the ABM to the path. 

simname = "RUN"; %define name of video for saving

timestep=6;%Indicate the amount of hours you want one timestep to take 
% Be careful not to increase timestep above 12 ({summed} probabilities will
% be higher than 1 in some cases, not doing anything anymore)

savePath = sprintf('output/ABM_run_variable/%s_t%dh/',simname,timestep);
if ~(exist(savePath,'dir'))
    mkdir(savePath); % create output directory to save results
end

[mySystem,cnst] = getSystemParams_variable([125,125],timestep);   %get system parameters
cnst.video = true;                             %show video on screen

%Adjust initial cell numbers (optional)               
mySystem.params.M1cellNo = 0; %starting number of M1 macrophages
mySystem.params.M2cellNo = 0; %starting number of M2 macrophages

%Additional modeling settings
mySystem.params.Fradius = 15; %Radius of the circle/ellips of the acinus ;
mySystem.params.Shape = 2;    %2 = it will be an ellips, 1 = it will be a circle
mySystem.params.Fellips = 15; %Radius of the second parameter defining ellips geometry (y-direction)
mySystem.params.Number = 1;   %Number of acini that you want to plot on the grid 

%run model 
[mySystemnew,TUcells,~,finalImage,TSC,Tm,Tm2,MutCells] = growTumor_variable(mySystem,cnst);

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
writeMyVideo(finalImage,"output/ABM_run_variable/vid_"+simname,4) %Create video for saving
saveas(gcf,savePath+"plot_"+simname+".png") %Save the video using the earlier provided name
