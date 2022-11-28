function [L, Fcells, Fprop,TUcells, TUprop,degd2,M1cells, M1prop, M2cells, M2prop] = CAF_go_grow_die(L, nh, Fcells, Fprop, Fpprol, Fpmig, Fpdeath, Frwalk, ChtaxMap,TUcells,TUprop,CFemt,ChtaxmapCF,CFprom,ChtaxMapS,CFmmp,ChtaxMapM1,CFmdiff,M1cells, M1prop, M2cells, M2prop,BM)
%Function in which CAFs perform their actions (migrate/proliferate/die). To incorporate directed migration, also
%the chemotaxis map for tumor cells is provided. Additionally, the actions of eliciting EMT in mutated cells, stimulating differentation of
%macrophages (from M1 to M2), creating a path in the stroma and promoting proliferation of tumor cells. The function returns an updated list of
%CAFs, macrophages, the stroma and tumor cells 

%CREATE MASKS FOR ADJACENT POSITIONS
m=getAdjacent(L,Fcells,nh);

%DETERMINE WHICH CELLS WILL PERFORM WHICH ACTIONS
P = m.randI2 <= Fpprol & ismember(m.S2,m.indxF); % GROW: proliferation (only if there is an empty spot in the neighborhood of the agent)
D = (m.randI2 <= (Fpdeath+Fpprol)) & (m.randI2 > Fpprol); % DIE: spontaneous death can happen any time
Mi = (m.randI2 <= (Fpdeath+Fpprol+Fpmig)) &  (m.randI2 > (Fpdeath+Fpprol)) & (ismember(m.S2,m.indxF)); % GO: migrate (only if there is an empty spot in the neighborhood of the agent)

De = P & (Fprop.Pcap(m.S2) == 0); % Also if the proliferation capacity is exhausted, the cell is added to the list of cells that will die this round
del = D | De; % cells to delete (due to spontaneous death and due to exhaustion of the proliferation probability)
act = find((P | Mi) & ~del); % indices to the cells that will perform action (migration or proliferation)

%ADD/REMOVE ACTING CELLS TO GRID
for iloop = 1:numel(act) % only for those that will do anything
    currID = act(iloop); % number within stack of currently acting cell
    ngh = m.S(:,m.S2(currID)); % cells neighborhood
    ngh2 = ngh(ngh>0); % erasing places that were previously occupied (cells cannot move to occupied spaces)
    indOL = find(~L(ngh2)); %selecting all free spots  
    chemo = ChtaxMap(ngh2(:)); % extract chemotaxis value (for tumor cells) at neighbors
    chemo = chemo(~L(ngh2)); % block occupied spots   
    if ~isempty(chemo) % use spot with highest chemo value
        chemo = chemo/max(chemo(:)); % normalize
        chemo = (1-Frwalk) * chemo + Frwalk * rand(size(chemo)); %Add the random effect for the migration action
        %[~,cid] = min(chemo); % lowest distance to tumor cells 
        [~,cid] = max(chemo); % greatest distance to tumor cells 
        indO = indOL(cid(1));   
        if ~isempty(indO) %if there is still a free spot
            L(ngh2(indO)) = true; %Set the selected grid space to occupied 
            if P(currID) % In case of proliferation
                Fcells = [Fcells uint32(ngh2(indO))]; % add new cell to stack
                Fprop.Pcap(m.S2(currID)) = Fprop.Pcap(m.S2(currID))-1; % decrease remaining prol. cap. of acting cell
                Fprop.Pcap = [Fprop.Pcap, Fprop.Pcap(m.S2(currID))]; % update property vector for Pmax (add the new cell)
            else % In case of migration
                L(Fcells(m.S2(currID))) = false; %Freeing the spot in which the cell was previously
                Fcells(m.S2(currID)) = uint32(ngh2(indO)); %Update the cells location in the list of CAFs 
            end
        end
    end
end

%CAFs eliciting EMT in tumor cells 
candidates = ChtaxmapCF(TUcells)<=3 & TUprop.EMT==0 & TUprop.mut>0; %Select tumor cells that are within reach and have already mutated 
numbers=rand(1,length(candidates)); %Generate random numbers to determine which cells will go through EMT 

if length(BM)>3
    if sum(candidates(:)) & BM(length(BM))<BM(2) %Tumor cells can only go through EMT after disease has become invasive (at least one basement membrane agent has been removed)
        Pdiff = candidates & (numbers <= CFemt); %Select the tumor cells that will go through EMT
        TUprop.EMT(Pdiff)=1; %Update the EMT property for those cells 
    end
end

%CAFs affecting macrophage differentiation (from M1 to M2)
m=getAdjacent(L,Fcells,nh); 
candidates = ChtaxMapM1(Fcells)<=5; %Only look for CAFs that are within a 5 gridspace distance of a macrophage 
if sum(candidates(:))
    Pdiff = candidates & (m.randI2 <= CFmdiff); %Select CAFs that will affect macrophage differentation based on the randomly generated number
    Pd2 = find(Pdiff==1);
    Pd=Fcells(Pd2); %select the CAFs from the list that will affect macrophage differentiation 
    N = 125;
    radius = 5;
    if ~isempty(Pd)
        for n = 1:numel(Pd) %Only for the CAFs that will affect macrophage differentation 
            %Select all gridspaces within a five gridspace radius of the acting CAF
            ngh = Pd(n)-radius:Pd(n)-1; 
            ngh = [ngh,Pd(n)+1:Pd(n)+radius];
            for nn = 1:radius
                ngh = [ngh,(Pd(n)-nn*N)-radius:(Pd(n)-nn*N)+radius];
                ngh = [ngh,(Pd(n)+nn*N)-radius:(Pd(n)+nn*N)+radius]; 
            end
            ngh2=ngh(ngh>0); 
            chemo = ChtaxMapM1(ngh2); %Determine distance to macrophages 
            chemo2 = chemo(chemo==0); %Select the spaces at which there is a macrophage (distance to it is 0)
            ngh3 = ngh2(chemo==0); 
            if ~isempty(chemo2) & sum(M2cells==ngh3(1))<1 %If there is a macrophage that it can affect and this is not a M2 macrophage
                id = find(M1cells==ngh3(1)); %select macrophage that will differentiate towards the M2 phenotype
                M1cells(id) = []; %Remove the M1 macrophage from the list
                M2prop.Kcap=[M2prop.Kcap,M1prop.Kcap(id)]; %Add the killing capacity property of the new M2 to the list
                M2prop.engaged=[M2prop.engaged,M1prop.engaged(id)]; %Add the engaged property of the new M2 to the list
                M1prop.Kcap(id) = [];      % remove Kmax property of the differentiated macrophage
                M1prop.engaged(id) = [];   % remove engagement property of the differentiated macrophage
                M2cells=[M2cells,ngh3(1)]; %Add the newly differentiated macrophage to the M2 list. 
            end
        end
    end
end

if ~isempty(del) % updating CAF death
    L(Fcells(m.S2(del))) = false;     % remove dead CAFs from grid
    Fcells(m.S2(del)) = [];            % remove dead CAFs from stack
    Fprop.Pcap(m.S2(del)) = [];      % remove Kmax property of dead CAFs
end

 %Allow for proliferation promotion by CAFs of tumor cells 
[TUcells, TUprop, Fcells, L] = CF_promote(TUcells, TUprop, Fcells, L,CFprom,nh,ChtaxMap);

%CAFs affecting permeability of the stroma (making space for invading cells)
degd2=[];
candidates2 = ChtaxMapS(Fcells) <=1; %Only select CAFs that are next to stroma
nrs = rand(1,length(candidates2)); %Generate random numbers for determination of what CAF will take action 
Pbms = (CFmmp <= nrs);  %Select acting CAFs
if sum(candidates2(:)) 
    Pd2 = find(Pbms==1); 
    Pd=Fcells(Pd2); %Select acting CAFs from the list 
    N = 125;
    radius = 6;
    if ~isempty(Pd)
        for n = 1:numel(Pd) %only for the CAFs that will act
            %determine what stromal agents are within the radius 
            ngh = Pd(n)-radius:Pd(n)-1; 
            ngh = [ngh,Pd(n)+1:Pd(n)+radius];
            for nn = 1:radius
                ngh = [ngh,(Pd(n)-nn*N)-radius:(Pd(n)-nn*N)+radius];
                ngh = [ngh,(Pd(n)+nn*N)-radius:(Pd(n)+nn*N)+radius]; 
            end
            %Do not allow for stroma removal at the boundaries of the grid
            rm1 = 1:125:(125*125);
            rm2 = 125:125:(125*125);
            degd2(ismember(degd2,rm1))=[];
            degd2(ismember(degd2,rm2))=[];
            degd2=[degd2,ngh(126<ngh<15499)];%select all remaining stromal spaces that are affected 
        end
        degd2 = unique(degd2);
    end
else
    degd2 = [];
end

end
