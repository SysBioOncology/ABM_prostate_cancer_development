function [L,Mcells,Mprop,TUcells,TUprop,Lt,killst] = mCellRound(L,Lt,Mcells,TUcells, Mprop,TUprop,MinfluxProb,MinfluxRate,Mkmax,Mspeed,Mpmig,Mrwalk,Mpkill,Mpdeath,MengagementDuration,ChtaxMap,nh,cell3,cell4,nh2)
%Round for the M1 macrophages to perform their actions; they can migrate, die and kill tumor cells. The output of the function is the updated list
%of macrophages and tumor cells. Additionally influxing macrohpages are added to the system during this round  

killst=0;   %Reset kill number of tumor cells
if MinfluxRate > 0 %Check if there even is an influxrate
if rand () <= MinfluxProb %Check if there is an influx of cells
    if L(253) ==0 %Check if the 'blood vessel entry' is free
        L(253) = true; %Add the influxing macrophage to the grid. 
        Mcells = [Mcells, 253]; %Add the new macrophage to the list of macrophages
        %Add the new macrophage's properties
        Mprop.Kcap = [Mprop.Kcap, repmat(Mkmax,1,1)]; 
        Mprop.engaged = [Mprop.engaged,zeros(1,1)];
    end
end
end

[Mcells,Mprop] = shuffleCells(Mcells,Mprop);    %shuffle cells
if numel(Mcells) > 0 %if there are any macrophages
for j = 1:(Mspeed-1)  %allow macrophages to move Mspeed-1 times per round
    L([Mcells,TUcells,cell3,cell4]) = true;  %ensure all cells are on the grid 
    [L, Mcells] = M_go2(Mcells, Mpmig, Mrwalk, ChtaxMap, L, nh2); %macrophage movement (with invasion of the acinus)
    [TUcells, TUprop, Mcells, Mprop, L, Lt,kills]  = M_kill (TUcells, TUprop, Mcells, Mprop, L, Lt, Mpkill, nh, ChtaxMap, MengagementDuration); %macrophage killing of tumor cells 
    killst=killst+sum(kills); %Keep track of the number of killed tumor cells by macrophages 
    Mprop.engaged(Mprop.engaged > 0) = Mprop.engaged(Mprop.engaged > 0)-1; %un-engage macrophage
end

%allow macrophages to move once more or die
[L, Mcells, Mprop] = M_go_die(Mcells, Mprop, Mpmig, Mpdeath, Mrwalk, ChtaxMap, L, nh);

end
   
end