function [L, TUcells, TUprop, TUdeath,Bcells,Lb] = TU_go_grow_die(L, nh, TUcells, TUprop, TUpprol, TUpmig, TUpdeath, TUps, TUrwalk, ChtaxMap,ChtaxMapC,width,ChtaxMapB,TUpkill,Bcells,Lb,TUthrshBM,For,BM)
%FUnction in which the luminal cells perform their actions: die, migrate, proliferate or idle. Additionally, it is the round during which tumor
%cells can break down the basement membrane (if they have gained enough mutations)

%CREATE MASKS FOR ADJACENT POSITIONS
m=getAdjacent(L,TUcells,nh);

%Create map for mutated tumor cells and accordingly changed pprol
factor = (1/1000);
MapAnt = TUprop.mut*factor+ TUpprol; %Increase proliferation probability for mutated cells 
MapAnt(TUprop.mut==0)=MapAnt(TUprop.mut==0)*0.2;
Mapdeath = ones(size(TUprop.mut))*TUpdeath; %Additionally alter the map for cell death of mutated luminal cells 
Mapdeath(TUprop.mut==0)=Mapdeath(TUprop.mut==0)*1;

%DETERMINE WHICH CELLS WILL PERFORM WHICH ACTIONS
Pres = m.randI2 <= MapAnt & ismember(m.S2,m.indxF); % GROW: proliferation (only if there is an empty spot in the neighborhood of the agent)
D = (m.randI2 <= (Mapdeath+MapAnt)) & (m.randI2 > MapAnt); % DIE: spontaneous death can happen any time
Mi = (m.randI2 <= (Mapdeath+TUpmig+MapAnt)) &  (m.randI2 > (Mapdeath+MapAnt)) & (ismember(m.S2,m.indxF)) &(TUprop.EMT>0); % GO: migrate (only if there is an empty spot in the neighborhood of the agent)
K = (m.randI2 <= (Mapdeath+TUpmig+MapAnt+TUpkill) & (m.randI2>(Mapdeath+MapAnt+TUpmig)) & (TUprop.mut>TUthrshBM)); %KILL: kill basement membrane agent 

Ps = Pres & rand(1,m.nC2) <= TUps & TUprop.isStem(m.S2) & ismember(m.S2,m.indxF); % symmetric division probability for stem cells 

De = Pres & (TUprop.Pcap(m.S2) == 0); % Also if the proliferation capacity is exhausted, the cell is added to the list of cells that will die this round
actt = Pres | Mi; % indices to the cells that will perform action (migration or proliferation)
del = D | De; % find dead / dying cells
act = find(actt & ~del); % live cells that will proliferate or migrate

%ADD/REMOVE ACTING CELLS TO GRID
for iloop = 1:numel(act) % only for those that will do anything
    currID = act(iloop); % number within stack of currently acting cell
    ngh = m.S(:,m.S2(currID)); % cells neighborhood
    %Remove neighbors that fall outside of the grid
    ngh2 = ngh(ngh>125);
    ngh2(ngh2>(125*125))=[];
    if length(BM)>3
        if BM(length(BM))==BM(2) %Check if the basement membrane has already been broken down (invasive disease)
            [x,y] = intersect(ngh2,For); %In that case, the invading cell will move outside of its original region (so out of the acinus)
            ngh2(y)=[];
        end
    end
    indOL = find(~L(ngh2)); %selecting all free spots  
    chemo = ChtaxMap(ngh2(:)); % extract chemotaxis value (of fibroblasts) at neighbors
    chemo = chemo(~L(ngh2)); % block occupied spots
    chemo2 = ChtaxMapC(ngh2(:)); % extract chemotaxis value (of basal cells) at neighbors
    chemo2 = chemo2(~L(ngh2)); %block occupied spots
    chemo3 = ChtaxMapB(ngh2(:)); % extract chemotaxis value (of basement membrane) at neighbors
    chemo3= chemo3(~L(ngh2)); %block occupied spots 
    if TUprop.mut(m.S2(currID)) == 0 %Only allowing non-mutated cells to proliferate within a certain distance from basal cell layer
        chemo = chemo(chemo2<=width & chemo3>1 & chemo3<3);
        indOL = indOL(chemo2<=width & chemo3>1 & chemo3<3);
    end
    if ~isempty(chemo) % if there is still a free spot
        chemo = (1-TUrwalk) * chemo + TUrwalk * rand(size(chemo)); %Add the random effect for the migration action
        [~,cid] = min(chemo); % lowest distance 
        indO = indOL(cid(1));   
        if ~isempty(indO) %if there is still a free spot
            L(ngh2(indO)) = true; %Set the selected grid space to occupied 
            if Pres(currID) % In case of proliferation
                newCell = uint32(ngh2(indO)); % find place for new cell
                TUcells = [TUcells, newCell(1)]; % add new cell to stack
                if Ps(currID) % symmetric division (new cell is a stem cell)
                    %Add properties of the new cell to the lists 
                    TUprop.isStem = [TUprop.isStem, true];
                    TUprop.Pcap = [TUprop.Pcap, TUprop.Pcap(m.S2(currID))]; 
                    TUprop.mut = [TUprop.mut, TUprop.mut(m.S2(currID))];
                    TUprop.EMT = [TUprop.EMT, TUprop.EMT(m.S2(currID))];
                else % new cell is not a stem cell
                    %Add properties of the new cell to the list
                    TUprop.isStem = [TUprop.isStem, false];
                    TUprop.Pcap = [TUprop.Pcap, TUprop.Pcap(m.S2(currID))-1]; 
                    TUprop.mut = [TUprop.mut, TUprop.mut(m.S2(currID))];
                    TUprop.EMT = [TUprop.EMT, TUprop.EMT(m.S2(currID))];
                     if ~TUprop.isStem(m.S2(currID)) & (TUprop.mut(m.S2(currID))>0)  % reduce proliferation capacity                     
                       TUprop.Pcap(m.S2(currID)) = TUprop.Pcap(m.S2(currID))-1;
                    end
                    if TUprop.mut(m.S2(currID)) & (TUprop.mut(m.S2(currID))>0) & (TUprop.Pcap(m.S2(currID))==1) > 0
                        TUprop.Pcap(length(TUprop.Pcap))=TUprop.Pcap(length(TUprop.Pcap))+1;
                    end
                end
            else % In case of migration
                L(ngh2(indO)) = false;
                [~,cid] = max(chemo); % lowest distance 
                indO = indOL(cid(1));
                L(ngh2(indO)) = true;
                L(TUcells(m.S2(currID))) = false; %Freeing the spot in which the cell was previously
                TUcells(m.S2(currID)) = uint32(ngh2(indO)); % update cell position in the list of tumor cells
            end
        end
    end
end

if sum(K)>0
    [Bcells, L,Lb] = TU_kill(Bcells, L,K,m,ChtaxMapB,Lb,TUcells); %Tumor cell kills basement membrane agent
end

if ~isempty(del) % remove dead tumor cells
    L(TUcells(m.S2(del))) = false;      % remove from grid
    TUcells(m.S2(del)) = [];            % remove from stack
    TUprop.isStem(m.S2(del)) = [];      % remove stemness property
    TUprop.Pcap(m.S2(del)) = [];        % remove Pmax property
    TUprop.mut(m.S2(del)) = [];         % remove mutation property
    TUprop.EMT(m.S2(del)) = [];         % remove EMT property
end 

TUdeath = sum(D);                        % number of dying tumor cells
end
