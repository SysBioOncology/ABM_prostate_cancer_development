function [L, Fcells, Fprop,CFcells,CFprop] = F_go_grow_die(L, nh, Fcells, Fprop, Fpprol, Fpmig, Fpdeath, Frwalk, ChtaxMap,CFcells,CFprop,Fdiff)
%Function in which fibroblasts perform their actions (migrate/proliferate/die). To incorporate directed migration, also
%the chemotaxis map for tumor cells is provided. However, as described in the paper, fibroblasts are mainly quiescent.
%So there will not be performed many actions. This can be altered in this script. The function returns the updated list of fibroblast and its properties
%It also returns the updated list of CAFs, as in this script, fibroblasts
%can also differentiate towards the CAF phenotype. 

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
    if ~isempty(chemo) 
        chemo = chemo/max(chemo(:)); % normalize
        chemo = (1-Frwalk) * chemo + Frwalk * rand(size(chemo)); %Add the random effect for the migration action
        [~,cid] = min(chemo); % lowest distance 
        indO = indOL(cid(1));   
        if ~isempty(indO) %if there is still a free spot
            L(ngh2(indO)) = true; %Set the selected grid space to occupied 
            if P(currID) % In case of proliferation
                Fcells = [Fcells uint32(ngh2(indO))]; % add new cell to stack
                Fprop.Pcap(m.S2(currID)) = Fprop.Pcap(m.S2(currID))-1; % decrease remaining prol. cap. of acting cell
                Fprop.Pcap = [Fprop.Pcap, Fprop.Pcap(m.S2(currID))]; % update property vector for Pmax (add the new cell)
            else % In case of migration
                L(Fcells(m.S2(currID))) = false; %Freeing the spot in which the cell was previously
                Fcells(m.S2(currID)) = uint32(ngh2(indO)); %Update the cells location in the list of fibroblasts 
            end
        end
    end
end

if ~isempty(del) % updating Fibroblast cell death
    L(Fcells(m.S2(del))) = false;     % remove dead Fibroblast from grid
    Fcells(m.S2(del)) = [];            % remove from stack
    Fprop.Pcap(m.S2(del)) = [];      % remove Pmax property
end
 
%Change fibroblasts into CAFs (if close enough to tumor cells)
candidates = ChtaxMap(Fcells)<=2; %Select fibroblasts that are close to mutated cells 

if sum(candidates(:)) %If there are fibroblasts that are close enough to mutated cells
    Pdiff = candidates & (m.randI2<Fdiff); %Draw random numbers and determine which cells will differentiate
    list = find(Pdiff==1);
    list=flip(list); %Go through the changing fibroblasts from end to front (so that you do not select the wrong cells due to changed indices)
    for k=1:length(list) %Go through the list of fibroblasts that will change
        j=list(k);
        CFprop.Pcap=[CFprop.Pcap,Fprop.Pcap(j)]; %Add the recently changed cell to the CAF property list
        Fprop.Pcap(j) = []; %Remove the property of the differentiated fibroblast from the list of fibroblasts
        CFcells=[CFcells,Fcells(j)]; %Add the recently changed cell to the CAF list
        Fcells(j) = []; %Remove the changed cell from the list of fibroblasts 
    end
end

end
