function mse = MSE(x,mySystem, cnst,  varNames, version, t, m)
%Calculate the Mean Squared Error for Particle Swarm Optimization
%Determines the MSE between data and simulation using a fitted parameter
%and returns this MSE. 

time = 0:4:cnst.nSteps*4; %Determine the time in the simulation to compare to time of measurement

%Adjust mySystem parameters for input parameters in x
%Loop to automatically change the parameters based on their name, should be
%adjusted a bit more before use. 
for i = 1:length(varNames)
    varName = varNames(i); %Should be: mySystem.params.varName 
    val = x(i);
end

%Manual adaptation of parameters to fit
TUpprol = x(1); %Change parameter 1
TUpmax = x(2); %Change parameter 2 


for j = 1:cnst.simnumber
    [~,TUcellNo,FcellNo,~] = growTumor(mySystem,cnst);                %run the ABM to determine simulation output 
    TUcellNo=TUcellNo(:,1);
    if version == "TU"  %LNCAP case
        CellnoRelative = 1+((TUcellNo-TUcellNo(1))/(TUcellNo(1)));    %calculate relative cell number of tumor cells 
    elseif version == "F" %Fibroblast case 
        CellnoRelative = 1+((FcellNo-FcellNo(1))/(FcellNo(1)));       %Calculate relative number of fibroblasts 
    else
        disp("Wrong version selected")
    end
    CellnoInt (:,j) = interp1(time,CellnoRelative,t);               %interpolate relative cell number to match time steps of data 
end

CellnoMean= mean(CellnoInt,2);                                   %calculate mean relative cell number


mse = mean((m - CellnoMean).^2);                               %calculate MSE between incucyte and model output

end