function [L, Mcells, Mprop] =  M_go_die(Mcells, Mprop, Mpmig, Mpdeath, Mrwalk, ChtaxMap, L, nh)
%At the end of their round, macrophages can migrate one more time or die(or do nothing). The function returns the updated list of macrophages and
%the updated grid. 
    
%CREATE MASKS FOR ADJACENT POSITIONS
m = getAdjacent(L,Mcells,nh);

%DETERMINE WHICH CELLS PERFORM ACTIONS
D = m.randI2 <= (Mpdeath); % DIE: spontaneous death can happen any time
Mi = (m.randI2 <= (Mpdeath+Mpmig)) &  (m.randI2 > (Mpdeath)) & (ismember(m.S2,m.indxF)); % GO: migrate migrate (only if there is an empty spot in the neighborhood of the agent)

del = D; % cells to delete
act = find(Mi & ~del); % indices to the cells that will perform action

%MOVE CELLS ON GRID
for iloop = 1:numel(act) % only for those that will do anything
    currID = act(iloop); % number within stack of currently acting cell
    ngh = m.S(:,m.S2(currID)); % cells neighborhood
    ngh2 = ngh(ngh>0); % erasing places that were previously occupied
    indOL = find(~L(ngh2)); %selecting all free spots  
    chemo = ChtaxMap(ngh2(:)); % extract chemotaxis value (of tumor cells) at neighbors
    chemo = chemo(~L(ngh2)); % block occupied spots   
    if ~isempty(chemo) % use spot with highest chemo value
        %Only if a macrophage is within reach of a mutated cell will it
        %migrate towards the tumor cell. Otherwise it will migrate randomly
        if min(chemo) <17
            Mrwalk = 0.8;
        else
            Mrwalk = 1;
        end
        chemo = chemo/max(chemo(:)); % normalize
        chemo = (1-Mrwalk) * chemo + Mrwalk * rand(size(chemo)); %Add the random effect for the migration action
        [~,cid] = min(chemo); % lowest distance 
        indO = indOL(cid(1));   
        if ~isempty(indO) %if there is still a free spot
            L(ngh2(indO)) = true; %Set the selected grid space to occupied 
            L(Mcells(m.S2(currID))) = false; %Freeing the spot in which the cell was previously
            Mcells(m.S2(currID)) = uint32(ngh2(indO)); %Update the cells location in the list of macrophages 
        end
    end
end

if ~isempty(del) % updating macrophage cell death
    L(Mcells(m.S2(del))) = false;    % remove macrophage cell from grid
    Mcells(m.S2(del)) = [];          % remove from stack
    Mprop.Kcap(m.S2(del)) = [];      % remove Kmax property
    Mprop.engaged(m.S2(del)) = [];   % remove engagement property
end
 
end
