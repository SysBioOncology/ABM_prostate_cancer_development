function [TUcells, TUprop, CFcells,L] = CF_promote(TUcells, TUprop, CFcells,L,CFprom,nh,ChtaxMap)
%This function allows for promotion of proliferation of tumor cells by CAFs
%It outputs the updated CAF and tumor cell list 

%Create masks for adjacent positions 
m = getAdjacent(L,TUcells,nh);

% pre-select tumor cells that are close enough (adjacent) to CAFs
candidates = ChtaxMap(CFcells)<=1;
if sum(candidates(:)) % if there are candidates
    % select cells that are going to promote
    K = candidates & (rand(1,length(CFcells))<CFprom);
    actK = find(K); % cell indices
    if ~isempty(actK) % if there is a cell that is going to promote
    targetIDs = int32(zeros(1,0)); % preallocate
    promoterIDs = int32(zeros(1,0)); % preallocate
    % start tumor cell promotion, same random order as before
    St = bsxfun(@plus,CFcells(actK),nh.aux(nh.Pms(:,randi(nh.nP,1,length(actK)))));
    % iterate through all CAFs and look at their neighborhood
    for jj = 1:size(St,2) 
        neighbPosit = St(randperm(length(nh.aux)),jj);
        instaprom = ismember(neighbPosit(:),TUcells(:));
        possibleTargets = neighbPosit(instaprom)'; % possible targets
        targetIDs=[targetIDs,int32(possibleTargets)]; %target cells to promote
        promoterIDs = [promoterIDs, CFcells(actK(jj))]; % add promoter ID to stack
    end

    % find indices to promoted cells and promoters. 
    auxPromTU = ismember(TUcells,targetIDs) & TUprop.mut>0.05; % which tumor cells are promoted
    act = find(auxPromTU==1);

    if sum(auxPromTU)>0 % if promoting happens, then check if promoted cells have space & let them proliferate
        for iloop=1:length(act) % for cells that can be promoted
            currID = act(iloop); % number within stack of currently acting cell
            ngh = m.S(:,m.S2(currID)); % cells neighborhood
            ngh2 = ngh(ngh>0);  % erasing places that were occupied
            indO = find(~L(ngh2),1,'first');  % selecting free spot
            if ~isempty(indO) % if there is still a free spot
                L(ngh2(indO)) = true; % add cell to grid
                newCell = uint32(ngh2(indO)); % find place for new cell
                TUcells = [TUcells, newCell(1)]; % add new cell to stack
                %Set properties of the new cell 
                TUprop.isStem = [TUprop.isStem, false]; 
                TUprop.Pcap = [TUprop.Pcap, TUprop.Pcap(m.S2(currID))-1];
                TUprop.mut = [TUprop.mut, TUprop.mut(m.S2(currID))];
                TUprop.EMT = [TUprop.EMT, TUprop.EMT(m.S2(currID))];
                if ~TUprop.isStem(m.S2(currID)) % reduce proliferation capacity
                    TUprop.Pcap(m.S2(currID)) = TUprop.Pcap(m.S2(currID))-1;
                end    
            end
        end
    end
    end 
end

end