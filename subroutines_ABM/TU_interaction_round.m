function [L, M1cells, M1prop, M2cells, M2prop,degd] = TU_interaction_round(L, nh, TUcells, TUpMdiff, TUthrshM, M1cells, M1prop, M2cells, M2prop, ChtaxMapM1,TUprop,TUpMMP,TUthrshMMP,degd,ChtaxMapS)
%Function that models the break down of the stroma by luminal cells and affected macrophage differentiation (from M1 macrophages to M2
%macrophages)

%CREATE MASKS FOR ADJACENT POSITIONS
m=getAdjacent(L,TUcells,nh);

%Affect macrophage differentiation 
candidates = ChtaxMapM1(TUcells)<=5; %Select tumor cells that are close enough to the macrophages
if sum(candidates(:))
    Pdiff = candidates & (m.randI2 <= TUpMdiff) & TUprop.mut>TUthrshM; %Select the tumor cells that are going to affect macrophage differentiation
    Pd2 = find(Pdiff==1);
    Pd=TUcells(Pd2); %Select the acting tumor cells
    N = 125;
    radius = 5;
    if ~isempty(Pd)
        for n = 1:numel(Pd) %Only for the agents that will actually act
            %Select the grid spaces that the tumor cell can reach 
            ngh = Pd(n)-radius:Pd(n)-1;
            ngh = [ngh,Pd(n)+1:Pd(n)+radius];
            for nn = 1:radius
                ngh = [ngh,(Pd(n)-nn*N)-radius:(Pd(n)-nn*N)+radius];
                ngh = [ngh,(Pd(n)+nn*N)-radius:(Pd(n)+nn*N)+radius]; 
            end
            ngh2=ngh(ngh>0);
            chemo = ChtaxMapM1(ngh2); %Determine distances to the macrophages 
            chemo2 = chemo(chemo==0); %Select the gridspaces that are occupied by a macrophage
            ngh3 = ngh2(chemo==0); 
            if ~isempty(chemo2) & sum(M2cells==ngh3(1))<1 %If there are macrophages that can react and not for M2 macrophages
                id = find(M1cells==ngh3(1)); %Select the differentiating macrophage
                M1cells(id) = []; %Remove the M1 macrophage from the stack of cells
                %Add the properties for the newly differentiated M2 macrophage
                M2prop.Kcap=[M2prop.Kcap,M1prop.Kcap(id)]; 
                M2prop.engaged=[M2prop.engaged,M1prop.engaged(id)];
                %Remove the properties for the differentiated M1 macrophage
                M1prop.Kcap(id) = [];      % remove Kmax property
                M1prop.engaged(id) = [];   % remove engagement property
                M2cells=[M2cells,ngh3(1)]; %Add the new M2 macrophage to the list
            end
        end
    end
end

%Degradation of surrounding stroma (to allow for invasion)
candidates2 = ChtaxMapS(TUcells) <=4; %Select tumor cells that are close enough to the stroma to affect it)
Pbms = (TUpMdiff <= m.randI2) & (m.randI2 <= (TUpMdiff + TUpMMP)) & TUprop.mut>TUthrshMMP & candidates2; %Determine which tumor cells are going to act
if sum(candidates2(:))
    Pd2 = find(Pbms==1);
    Pd=TUcells(Pd2); %Select the agents that are actually going to act 
    N = 125;
    radius = 6;
    if ~isempty(Pd) %If there are tumor cells that are going to affect the stroma
        for n = 1:numel(Pd) %Only for those that will act 
            %Select the area of stroma that could be affected by the acting tumor cells 
            ngh = Pd(n)-radius:Pd(n)-1;
            ngh = [ngh,Pd(n)+1:Pd(n)+radius];
            for nn = 1:radius
                ngh = [ngh,(Pd(n)-nn*N)-radius:(Pd(n)-nn*N)+radius];
                ngh = [ngh,(Pd(n)+nn*N)-radius:(Pd(n)+nn*N)+radius]; 
            end
            rm1 = 1:125:(125*125); %Part hard-coded (should be changed if you change the grid size)
            rm2 = 125:125:(125*125);
            degd(ismember(degd,rm1))=[]; %Remove the boundaries of the grid, as tumor cells should never be allowed to move here
            degd(ismember(degd,rm2))=[];
            degd=[degd,ngh(126<ngh<15499)]; %Remove all stromal agents that are affected by the tumor cells 
        end
        degd = unique(degd);
    end
end
end