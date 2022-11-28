 function [L, Mcells] = M_go2(Mcells, Mpmig, Mrwalk, ChtaxMap, L, nh2)
 %Migration of macrohpages (they can invade the acinus, and have an
 %enlarged neighborhood towards which they can migrate (twice the Moore
 %neighborhood). The output of this function is the updated macrophage list
 %and the updated grid. 

%CREATE MASKS FOR ADJACENT POSITIONS
m = getAdjacent2(L,Mcells,nh2); 

%DETERMINE WHICH CELLS WILL MIGRATE
Mi = m.randI <= Mpmig; % indices to cells that will migrate
act = find(Mi); % indices to the cells that will perform action

%MOVE CELLS ON GRID
for iloop = 1:numel(act) % only for those that will do anything
    currID = act(iloop); % number within stack of currently acting cell
    ngh = m.S(:,m.indxF(currID)); % cells neighborhood
    ngh2 = ngh(ngh>0); % erasing places that were previously occupied (cells cannot move to occupied spaces)
    indOL = find(~L(ngh2)); %selecting all free spots  
    chemo = ChtaxMap(ngh2(:)); % extract chemotaxis (for tumor cells) value at neighbors
    chemo = chemo(~L(ngh2)); % block occupied spots   
    if ~isempty(chemo) 
        %If within reach of a tumor cell, the macrophage moves towars the
        %tumor cell, otherwise it moves randomly
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
            L(Mcells(m.indxF(currID))) = false; %Freeing the spot in which the cell was previously
            Mcells(m.indxF(currID)) = uint32(ngh2(indO)); %Update the cells location in the list of macrophages
        end
    end
end

 
end