function [L, Ccells, Cprop] = C_go_grow_die(L, nh, Ccells, Cprop, Cpprol, Cpmig, Cpdeath, Crwalk, ChtaxMap,ChtaxMapB)
%Function in which basal cells perform their actions
%(migrate/proliferate/die). To incorporate directed migration, also
%chemotaxis maps are provided as input (of the fibroblasts and the basement
%membrane). The output is the updated list of basal cells 

%CREATE MASKS FOR ADJACENT POSITIONS
m=getAdjacent(L,Ccells,nh);

%DETERMINE WHICH CELLS WILL PERFORM WHICH ACTIONS
P = m.randI2 <= Cpprol & ismember(m.S2,m.indxF); % GROW: proliferation (only if there is an empty spot in the neighborhood of the agent
D = (m.randI2 <= (Cpdeath+Cpprol)) & (m.randI2 > Cpprol); % DIE: spontaneous death can happen any time
Mi = (m.randI2 <= (Cpdeath+Cpprol+Cpmig)) &  (m.randI2 > (Cpdeath+Cpprol)) & (ismember(m.S2,m.indxF)); % GO: migrate (only if there is an empty spot in the neighborhood of the agent)
%Agents will perform an action if the randomly generated number (from
%randI2 is in the defined range above (using the prescribed probability)

De = P & (Cprop.Pcap(m.S2) == 0); % Also if the proliferation capacity is exhausted, the cell is added to the list of cells that will die this round
del = D | De; % cells to delete (due to spontaneous death and due to exhaustion of the proliferation probability)
act = find((P | Mi) & ~del); % indices to the cells that will perform action (migration or proliferation)

%ADD/REMOVE ACTING CELLS TO GRID
for iloop = 1:numel(act) % only for those that will do anything
    currID = act(iloop); % number within stack of currently acting cell
    ngh = m.S(:,m.S2(currID)); % determine cells in its neighborhood
    ngh2 = ngh(ngh>0); % erasing places that were previously occupied (cells cannot move to occupied spaces)
    indOL = find(~L(ngh2)); %selecting all free spots  
    chemo = ChtaxMap(ngh2(:)); % extract chemotaxis value at neighbors (for fibroblasts)
    chemo = chemo(~L(ngh2)); % block occupied spots
    chemo2 = ChtaxMapB(ngh2(:)); % extract chemotaxis value at neighbors (for basement membrane)
    chemo2 = chemo2(~L(ngh2)); %block occupied spots 
    chemo = chemo(chemo2==1); %Only select spots that are next to a basement membrane agent(basal cells should be attached to basement membrane)
    indOL = indOL(chemo2==1);
    if ~isempty(chemo) % use spot with highest chemo value
        chemo = chemo/max(chemo(:)); % normalize
        chemo = (1-Crwalk) * chemo + Crwalk * rand(size(chemo)); %Add the random effect for the migration action
        [~,cid] = min(chemo); % lowest distance to the agent (fibroblast in this case)
        indO = indOL(cid(1));   
        if ~isempty(indO) %if there is still a free spot that meets all requirements 
            L(ngh2(indO)) = true; % Set the selected grid space to occupied 
            if P(currID) % In case of proliferation
                Ccells = [Ccells uint32(ngh2(indO))]; % add new cell to stack
                Cprop.Pcap(m.S2(currID)) = Cprop.Pcap(m.S2(currID))-1; % decrease remaining prol. cap. of acting cell 
                Cprop.Pcap = [Cprop.Pcap, Cprop.Pcap(m.S2(currID))]; % update property vector for Pmax (add the new cell)
            else % In case of migration
                L(Ccells(m.S2(currID))) = false; %Freeing the spot in which the cell was previously
                Ccells(m.S2(currID)) = uint32(ngh2(indO)); %Update the cells location in the list of basal cells 
            end
        end
    end
end

if ~isempty(del) % Performing basal cell death
    L(Ccells(m.S2(del))) = false;     % remove dead basal cells from grid
    Ccells(m.S2(del)) = [];            % remove from stack
    Cprop.Pcap(m.S2(del)) = [];      % remove Pmax property of dead cells
end
 
end