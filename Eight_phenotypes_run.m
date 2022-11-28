%Script to run multiple simulations for all eight phenotypes and plot the
%average behavior over time. 

CFpromm = [0.7 0.7 0.7 0.01 0.01 0.01 0.01 0.7]; %List of CAF promotion for the eight different phenotypes 
ThrsMm = [20 20 1 20 1 20 1 1]; %List of the mutation threshold for M1 to M2 differentiation for the eight different phenotypes 
TUpmutt = [0.000175 0.0000075 0.000175 0.000175 0.000175 0.0000075 0.0000075 0.0000075]; %List of mutation probability for the eight different phenotypes 
M1pmig =[0.01 0.01 0.51 0.01 0.51 0.01 0.51 0.51]; %List of M1 migration probability for the eight different phenotypes 

colors = ["black", "blue", "green", "yellow", "red", "magenta", "cyan","grey"];
additional = ["a","b","c","d","e","f","g","h"];
xx = 25; %Number of simulations per phenotype

%Add all the relevant folders to the path
addpath('./functions/') 
addpath('./subroutines_ABM/')


for i= 1:length(CFpromm)
    
      %Create lists to store variables for plotting
      tumcells=[];
      mutcells =[];
      timeinv = [];
      
    dim = 0;
    while dim(1)<xx
        %We only want to save the simulations in which disease progresses
        %towards cancer, which is why we use a while loop. 

        [mySystem,cnst] = getSystemParams_long([125,125]);   %get system parameters
        
        %Overwrite the parameters that change for the simulated phenotype
        %(select the value from the list) 
        mySystem.params.CFprom=CFpromm(i);
        mySystem.params.TUthrshM=mySystem.params.TUantig*ThrsMm(i);
        mySystem.params.TUpmut=TUpmutt(i);
        mySystem.params.M1pmig=M1pmig(i);
        
        cnst.video = false; %Turn off video, since it really slows the simulations if you want to run it many times
        mySystem.params.Fradius = 15; %Set radius of prostatic acinus 

        %run model 
        [mySystem,TUcellNo,FcellNo,finalImage,TSC,Tm,Tm2,MutCells,BM] = growTumor(mySystem,cnst);

        %Save relative tumor number and time of invasion 
        TUnoRelative = 1+((TUcellNo-TUcellNo(1))/(TUcellNo(1))); %calculate relative cell number 
        yy = find(BM==BM(2)-1); %See if there is a time during the past simulation in which basement membrane was broken down (invasion = cancer)
        
        %only if the simulated case developed to cancer, the information
        %will be saved
        if ~isempty(yy)
                tumcells = [tumcells ; TUnoRelative(:,1)']; %Save and store relative tumor cell numbers over time
                mutcells = [mutcells ; Tm2'];    %Save the absolute number of mutated cells
                timeinv = [timeinv; yy(1)];      %Save the timepoint at which the first part of basement membrane was broken down (the point at which it became cancer)
                %Additional parameters that you want to track over time can
                %be saved in this loop as well 
        end
        dim = size(tumcells); %Check the dimension; see how many simulations that resulted in cancer were already conducted for this phenotype 
    end
    sz=size(tumcells);
    time = linspace(0,cnst.nSteps*0.5,cnst.nSteps+1); %Create a time variable for plotting 
    
    %Plot some results 
    if sz(1)>1
        %Plot the relative nr of tumorcells over time for each phenotype 
        figure(1), hold on %hold on to plot all eight phenotypes in graph
        ploterr(time,tumcells,'Color',colors(i),'Error','sem') %this function plots the average and as a shade the standard deviation or standard error of the mean (user can choose)
        xlabel('time (days)'); ylabel('Relative amount of luminal cells'); title('Standard Error');
        
        %Plot the absolute nr of mutated cells for each phenotype
        figure(2), hold on
        ploterr(time,mutcells,'Color',colors(i),'Error','sem') %Plot the average behavior over time and shaded the std or standard error of the mean 
        xlabel('time (days)'); ylabel('Absolute amount of mutated cells');
    end
    
    %Create legends linking the phenotypes to the graphs
    figure(1), hold on
    legend('std','CFhigh M anti Pmut high','std','CFhigh Manti Pmut low','std','CFhigh M prom Pmut high','std','CFlow  M anti Pmut high','std','CFlow M prom Pmut high','std','CFlow Manti Pmut low','std','CFlow M prom Pmut low','std','CFhigh M prom Pmut low');
    
    figure(2), hold on
    legend('std','CFhigh M anti Pmut high','std','CFhigh Manti Pmut low','std','CFhigh M prom Pmut high','std','CFlow  M anti Pmut high','std','CFlow M prom Pmut high','std','CFlow Manti Pmut low','std','CFlow M prom Pmut low','std','CFhigh M prom Pmut low');
    
    save(colors(i),'mutcells') %Save the number of mutated cells for each timestep for all iterations for all phenotypes
    save(additional(i),'timeinv') %Save the times at which each simulation became invasive 
end
      