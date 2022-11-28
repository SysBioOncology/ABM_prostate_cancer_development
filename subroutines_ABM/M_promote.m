function [TUcells, TUprop, Mcells, Mprop, L, Lt] = M_promote(TUcells, TUprop, Mcells, Mprop, L, Lt,Mpprom,nh,ChtaxMap,engagementDuration)
%Function that lets macrophages promote tumor cell proliferation. It returns updated lists of macrophages, tumor cells and an updated grid. 

%CREATE MASKS FOR ADJACENT POSITIONS
m = getAdjacent(L,TUcells,nh);

% pre-select macrophages that may be close enough to tumor cells (adjacent)
candidates = ChtaxMap(Mcells)<=1;
if sum(candidates(:)) % if there are candidates
    % select cells that are going to promote
    K = candidates & (Mprop.engaged==0) & (Mprop.Kcap>0) & (rand(1,length(Mcells))<Mpprom);
    actK = find(K); % cell indices
    if ~isempty(actK) % if there is a cell that is going to promote
    targetIDs = int32(zeros(1,0)); % preallocate
    promoterIDs = int32(zeros(1,0)); % preallocate
    % start tumor cell promotion, same random order as before
    St = bsxfun(@plus,Mcells(actK),nh.aux(nh.Pms(:,randi(nh.nP,1,length(actK)))));
    % iterate through all macrophages and look at their neighborhood
    for jj = 1:size(St,2) 
        neighbPosit = St(randperm(length(nh.aux)),jj);
        instaprom = ismember(neighbPosit(:),TUcells(:));
        possibleTargets = neighbPosit(instaprom)'; % possible targets
        targetIDs=[targetIDs,int32(possibleTargets)]; %target cells to promote
        promoterIDs = [promoterIDs, Mcells(actK(jj))]; % add promoter ID to stack
    end

    % find indices to promoted cells and promoters. 
    auxPromTU = ismember(TUcells,targetIDs) & TUprop.mut>0.05; % which tumor cells are promoted
    auxPromM = ismember(Mcells,promoterIDs); % which macrophages are going to promote the tumor cells 
    act = find(auxPromTU==1);

    if sum(auxPromTU)>0 % if promoting happens, then check if promoted cells have space & let them proliferate
        for iloop=1:length(act) % for cells that can be promoted
            currID = act(iloop); % number within stack of currently acting cell
            ngh = m.S(:,m.S2(currID)); % cells neighborhood
            ngh2 = ngh(ngh>126);  % erasing places that were occupied
            ngh2(ngh2>(125*123))=[]; %Do not allow for proliferation into boundaries or outside of the grid
            indO = find(~L(ngh2),1,'first');  % selecting free spot
            if ~isempty(indO) % if there is still a free spot
                L(ngh2(indO)) = true; %Set the selected grid space to occupied 
                newCell = uint32(ngh2(indO)); % find place for new cell
                TUcells = [TUcells, newCell(1)]; % add new cell to stack
                %Add properties for the new cell 
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
    Mprop.engaged(auxPromM) = engagementDuration; % promoters are engaged
    L(TUcells) = true;  % Ensure all tumor cells are on the grid 
    Lt(TUcells) = true;  % And also on the tumor grid
    end
end


end