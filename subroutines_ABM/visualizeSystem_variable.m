function visualizeSystem_variable(TUprop,mySystem,all)
%Function that visualizes the running simulation (it shows all agents and
%the current time step - for a variable time step)

%Load all agents and properties 
cellfun(@(f) evalin('caller',[f ' = TUprop.' f ';']), fieldnames(TUprop));
cellfun(@(f) evalin('caller',[f ' = mySystem.params.' f ';']), fieldnames(mySystem.params));
cellfun(@(f) evalin('caller',[f ' = mySystem.grid.' f ';']), fieldnames(mySystem.grid));
cellfun(@(f) evalin('caller',[f ' = mySystem.TU.' f ';']), fieldnames(mySystem.TU));
cellfun(@(f) evalin('caller',[f ' = mySystem.F.' f ';']), fieldnames(mySystem.F));
cellfun(@(f) evalin('caller',[f ' = mySystem.B.' f ';']), fieldnames(mySystem.B));
cellfun(@(f) evalin('caller',[f ' = mySystem.C.' f ';']), fieldnames(mySystem.C));
cellfun(@(f) evalin('caller',[f ' = mySystem.M1.' f ';']), fieldnames(mySystem.M1));
cellfun(@(f) evalin('caller',[f ' = mySystem.M2.' f ';']), fieldnames(mySystem.M2));
cellfun(@(f) evalin('caller',[f ' = mySystem.CF.' f ';']), fieldnames(mySystem.CF));

%RGB color coding 
F = ones(N,M,3);

%Plot stroma (define colors R-G-B)
F(all) = 255/255;
F(all+(N*M)) = 247/255;
F(all+(N*M)*2) = 171/255;

%PLOT NORMAL TUMOR CELLS
TUcellnormal = TUcells(~isStem);
F(TUcellnormal) = 120/255;
F(TUcellnormal+(N*M)) = 30/255;
F(TUcellnormal+(N*M)*2) = 30/255;

%PLOT MUTATED CELLS 
mut2 = mut>0;
TUcellmut = TUcells(mut2);
F(TUcellmut) = 218/255;
F(TUcellmut + (N*M)) = 23/255;
F(TUcellmut + (N*M)*2)= 23/255;

%PLOT MUTATED CELLS THAT WENT THROUGH EMT
EMT2 = EMT>0;
TUcellEMT = TUcells(EMT2);
F(TUcellEMT) = 226;
F(TUcellEMT + (N*M)) = 88/255;
F(TUcellEMT + (N*M)*2) = 62/255;

%PLOT STEM CELLS 
TUcellstem = TUcells(isStem);
F(TUcellstem) = 255/255;
F(TUcellstem + (N*M)) = 151/255;
F(TUcellstem + (N*M)*2)=98/255;

%PLOT MUTATED STEM CELLS 
TUcellmut2 = TUcells(mut2 & isStem);
F(TUcellmut2) = 230/255;
F(TUcellmut2 + (N*M)) = 165/255;
F(TUcellmut2 + (N*M)*2)=100/255;

%PLOT FIBROBLASTS
F(Fcells) = 61/255;
F(Fcells+(N*M)) = 129/255;
F(Fcells+(N*M)*2) = 14/255;

%CAFs
F(CFcells) = 102/255;
F(CFcells+(N*M)) = 204/255;
F(CFcells+(N*M)*2) = 0/255;

%PLOT COLUMNAR CELLS
F(Ccells) = 130/255;
F(Ccells+(N*M)) = 130/255;
F(Ccells+(N*M)*2) = 130/255;

%PLOT BASEMENT MEMBRANE
F(Bcells) = 77/255;
F(Bcells+(N*M)) = 77/255;
F(Bcells+(N*M)*2) = 77/255;

%PLOT M1
F(M1cells) = 15/255;
F(M1cells+(N*M)) = 109/255;
F(M1cells+(N*M)*2) = 214/255;

%PLOT M2
F(M2cells) = 50/255;
F(M2cells+(N*M)) = 192/255;
F(M2cells+(N*M)*2) = 255/255;

%PLOT FIGURE PROPERTIES
imshow(F,'InitialMagnification',600)
beginSc=100;
text(beginSc,beginSc+15,[num2str(round(mySystem.grid.StepsDone*mySystem.params.timestep)),' hours'],...
    'Color','k','FontWeight','bold','FontSize',18,'VerticalAlignment','top')
end


