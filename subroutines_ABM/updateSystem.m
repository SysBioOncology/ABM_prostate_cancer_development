function mySystem = updateSystem(mySystem,TUcells,Fcells,Bcells,Ccells,M1cells,M2cells,CFcells,TUprop,M1prop,M2prop,Fprop,Cprop,CFprop,L,Lt,Ls,i)
%Function that updates the system after each iteration. It returns an updated version of the system (all cell lists, properties, parameters and
%the grid)

%Copy variables back to mySystem
mySystem.TU.TUcells = TUcells;
mySystem.TU.TUprop.isStem = TUprop.isStem;
mySystem.TU.TUprop.Pcap = TUprop.Pcap;
mySystem.TU.TUprop.mut = TUprop.mut;
mySystem.TU.TUprop.EMT = TUprop.EMT;

mySystem.F.Fcells = Fcells;
mySystem.F.Fprop.Pcap = Fprop.Pcap;

mySystem.CF.CFcells=CFcells;
mySystem.CF.CFprop.Pcap=CFprop.Pcap;

mySystem.B.Bcells = Bcells;

mySystem.C.Ccells = Ccells;
mySystem.C.Cprop.Pcap = Cprop.Pcap;

mySystem.M1.M1cells = M1cells;
mySystem.M1.M1prop.Kcap = M1prop.Kcap;
mySystem.M1.M1prop.engaged = M1prop.engaged;

mySystem.M2.M2cells = M2cells;
mySystem.M2.M2prop.Kcap = M2prop.Kcap;
mySystem.M2.M2prop.engaged = M2prop.engaged;

%Update the grid 
mySystem.grid.L = L;
mySystem.grid.Lt = Lt;
mySystem.grid.Ls = Ls;

%Update the amount of iterations that was already conducted 
if isfield(mySystem.grid, 'StepsDone')
    mySystem.grid.StepsDone = mySystem.grid.StepsDone + 1;
else 
    mySystem.grid.StepsDone = i;
end

end