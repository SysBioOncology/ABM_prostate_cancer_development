%Script to run multiple simulations for all eight phenotypes and plot the
%average behavior over time. 

CFpromm = [0.7 0.7 0.7 0.01 0.01 0.01 0.01 0.7]; %List of CAF promotion for the eight different phenotypes 
ThrsMm = [20 20 1 20 1 20 1 1]; %List of the mutation threshold for M1 to M2 differentiation for the eight different phenotypes 
TUpmutt = [0.000175 0.0000075 0.000175 0.000175 0.000175 0.0000075 0.0000075 0.0000075]; %List of mutation probability for the eight different phenotypes 
M1pmig =[0.01 0.01 0.51 0.01 0.51 0.01 0.51 0.51]; %List of M1 migration probability for the eight different phenotypes 

colors = ["black", "blue", "green", "yellow", "red", "magenta", "cyan","grey"];
RGBad = [[0.5 0.5 0.5]]; %There are just eight default colors in Matlab, one of which is white. Here you can specify additonal colors using their RGB code
additional = ["a","b","c","d","e","f","g","h"]; %Vector to store the time at which the saved simulation became invasive
xx = 2;%10; %Number of simulations per phenotype (the higher, the more accurate the prediction will get, due to the stochasticity of the simulations).  

%Add all the relevant folders to the path
addpath('./functions/') 
addpath('./subroutines_ABM/')

savePath = 'output/Eight_phenotypes_run/';
if ~(exist(savePath,'dir'))
    mkdir(savePath); % create output directory to save results
    mkdir(strcat(savePath,'timeinvasive/'))
    mkdir(strcat(savePath,'mutcells/'))
end

for i= 1%:length(CFpromm)
    
    fprintf('Start phenotype %d\n',i)
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
        mySystem.params.TUcellNo = 0;                
        mySystem.params.M1cellNo = 0;
        mySystem.params.M2cellNo = 0;
        mySystem.params.Shape = 2; %2 = it will be an ellips, 1 = it will be a circle
        mySystem.params.Fellips = 15;
        mySystem.params.Number = 1;%3;
        

        %run model 
        [mySystem,TUcellNo,FcellNo,finalImage,TSC,Tm,Tm2,MutCells,BM] = growTumor(mySystem,cnst);

        %Save relative tumor number and time of invasion 
        TUnoRelative = 1+((TUcellNo-TUcellNo(1))/(TUcellNo(1))); %calculate relative cell number 
        yy = find(BM==BM(2)-1); %See if there is a time during the past simulation in which basement membrane was broken down (invasion = cancer)
        check = TUnoRelative(400:end,1);
        
        %only if the simulated case developed to cancer, the information
        %will be saved
        if ~isempty(yy) %& sum(check)>300
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
        if i<8
            figure(1), hold on %hold on to plot all eight phenotypes in graph
            ploterr(time,tumcells,'Color',colors(i),'Error','sem') %this function plots the average and as a shade the standard deviation or standard error of the mean (user can choose)
            xlabel('time (days)'); ylabel('Relative amount of luminal cells'); title('Tumor cells over time');
        
            %Plot the absolute nr of mutated cells for each phenotype
            figure(2), hold on
            ploterr(time,mutcells,'Color',colors(i),'Error','sem') %Plot the average behavior over time and shaded the std or standard error of the mean 
            xlabel('time (days)'); ylabel('Absolute amount of mutated cells'); title('Mutated tumor cells over time');
        else
            figure(1), hold on %hold on to plot all eight phenotypes in graph
            ploterr(time,tumcells,'Color',RGBad(i-7,:),'Error','sem') %this function plots the average and as a shade the standard deviation or standard error of the mean (user can choose)
            xlabel('time (days)'); ylabel('Relative amount of luminal cells'); title('Tumor cells over time');
            legend(["std","CFhigh M anti Pmut high","std","CFhigh Manti Pmut low","std","CFhigh M prom Pmut high","std","CFlow  M anti Pmut high","std","CFlow M prom Pmut high","std","CFlow Manti Pmut low","std","CFlow M prom Pmut low","std","CFhigh M prom Pmut low"],'Location','eastoutside');

            figure(2), hold on
            ploterr(time,mutcells,'Color',RGBad(i-7,:),'Error','sem') %Plot the average behavior over time and shaded the std or standard error of the mean 
            xlabel('time (days)'); ylabel('Absolute amount of mutated cells'); title('Mutated tumor cells over time');
            legend(["std","CFhigh M anti Pmut high","std","CFhigh Manti Pmut low","std","CFhigh M prom Pmut high","std","CFlow  M anti Pmut high","std","CFlow M prom Pmut high","std","CFlow Manti Pmut low","std","CFlow M prom Pmut low","std","CFhigh M prom Pmut low"],'Location','eastoutside');
        end
    end
    
    %Create legends linking the phenotypes to the graphs
    save(savePath+"mutcells/"+colors(i),'mutcells') %Save the number of mutated cells for each timestep for all iterations for all phenotypes
    save(savePath+"timeinvasive/"+additional(i),'timeinv') %Save the times at which each simulation became invasive 
end

saveas(figure(1),savePath+"tumcells.png")
saveas(figure(2),savePath+"mutcells.png")
