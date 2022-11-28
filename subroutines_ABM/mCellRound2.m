function [L,Mcells,Mprop,TUcells,TUprop,Lt] = mCellRound2(L,Lt,Mcells,TUcells, Mprop,TUprop,MinfluxProb,MinfluxRate,Mkmax,Mspeed,Mpmig,Mrwalk,Mpkill,Mpdeath,MengagementDuration,ChtaxMap,nh2,cell3,cell4,nh,M2TUadd,ChtaxMapM2,M2emt,BM)
%Round for the M2 macrophages to perform their actions; they can migrate, die, kill and promote tumor cell proliferation or EMT. The output of the function is the updated list
%of macrophages and tumor cells. Additionally influxing macrohpages are added to the system during this round  

%In the paper, the influxrate was set to 0 (M2 macrophage do not enter the simulation by influx, instead they result from differentiation of M1
%macrophages)
if MinfluxRate > 0            %Check if there even is an influxrate
if rand () <= MinfluxProb     %Check if there is an influx of cells
    if L(253) ==0 %Check if the 'blood vessel entry' is free
        L(253) = true; %Add the influxing macrophage to the grid
        Mcells = [Mcells, 253]; %Add the new macrophage to the list of macrophages
        %Add the new macrophage's properties
        Mprop.Kcap = [Mprop.Kcap, repmat(Mkmax,1,1)];
        Mprop.engaged = [Mprop.engaged,zeros(1,1)]; 
    end
end
end

[Mcells,Mprop] = shuffleCells(Mcells,Mprop);    %shuffle cells
if numel(Mcells) > 0  %if there are any macrophages
for j = 1:(Mspeed-1)  %allow macrophages to move Mspeed-1 times per round
    L([Mcells,TUcells,cell3,cell4]) = true;  %ensure all cells are on the grid
    [L, Mcells] = M_go2(Mcells, Mpmig, Mrwalk, ChtaxMap, L, nh2);  %macrophage movement (also allowing invasion of the acinus)
    [TUcells, TUprop, Mcells, Mprop, L, Lt,k]  = M_kill (TUcells, TUprop, Mcells, Mprop, L, Lt, Mpkill, nh, ChtaxMap, MengagementDuration); %macrophage killing of tumor cells 
    Mprop.engaged(Mprop.engaged > 0) = Mprop.engaged(Mprop.engaged > 0)-1; %un-engage macrophage
    [TUcells, TUprop, Mcells, Mprop, L, Lt] = M_promote(TUcells, TUprop, Mcells, Mprop, L, Lt, M2TUadd, nh, ChtaxMap, MengagementDuration); %Promote tumor cell proliferation 
end

%allow macrophages to move once more or die
[L, Mcells, Mprop] = M_go_die(Mcells, Mprop, Mpmig, Mpdeath, Mrwalk, ChtaxMap, L, nh);

%EMT if tumor cells are close enough
candidates = ChtaxMapM2(TUcells)<=3 & TUprop.EMT==0 & TUprop.mut>0; %Select tumor cells that are within reach and have already mutated 
numbers=rand(1,length(candidates)); %Generate random numbers to determine which cells will go through EMT

if length(BM)>3
    if sum(candidates(:)) & BM(length(BM))<BM(2) %Tumor cells can only go through EMT after disease has become invasive (at least one basement membrane agent has been removed)
        Pdiff = candidates & (numbers <= M2emt); %Select the tumor cells that will go through EMT
        TUprop.EMT(Pdiff)=1; %Update the EMT property for those cells 
    end
end

end
    
end