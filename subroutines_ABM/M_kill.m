function [TUcells, TUprop, Mcells, Mprop, L, Lt,auxKillTU] = M_kill(TUcells, TUprop, Mcells, Mprop, L, Lt,Mpkill,nh,ChtaxMap,engagementDuration)
%Function that performs the killing of tumor cells by the macrophages. It returns the updated lists of tumor cells and macrophages and the updated grid. 

% pre-select macrophages that are close enough (adjacent) to tumor cells
candidates = ChtaxMap(Mcells)<=1;
if sum(candidates(:)) % if there are candidates
    % select cells that are going to kill
    K = candidates & (Mprop.engaged==0) & (Mprop.Kcap>0) & (rand(1,length(Mcells))<Mpkill);
    actK = find(K); % cell indices
    if ~isempty(actK) % if there is a cell that is going to kill
    targetIDs = int32(zeros(1,0)); % preallocate
    killerIDs = int32(zeros(1,0)); % preallocate
    % start tumor cell killing, same random order as before
    St = bsxfun(@plus,Mcells(actK),nh.aux(nh.Pms(:,randi(nh.nP,1,length(actK)))));
    % iterate through all macrophages and look at their neighborhood
    for jj = 1:size(St,2) 
        neighbPosit = St(randperm(length(nh.aux)),jj);
        instakill = ismember(neighbPosit(:),TUcells(:));
        % if the cell encounters another cell to kill
        if sum(instakill) > 0 
            % if more than 1 possible targets then use the first one
            possibleTargets = neighbPosit(instakill); % possible targets
            killwhat = int32(possibleTargets(1)); % kill only the first candidate   
            targetIDs = [targetIDs, killwhat]; % add target ID to stack
            killerIDs = [killerIDs, Mcells(actK(jj))]; % add killer ID to stack
        end
    end

    % find indices to killed cell and killer cell. If the unlikely case
    % happens that one tumor cell is simultaneously killed by two macrophages,
    % then both will be exhausted
    auxKillTU = ismember(TUcells,targetIDs) & TUprop.mut>0; % which tumor cells are killed
    auxKillM = ismember(Mcells,killerIDs); % which macrophages perform the killing action 

    if sum(auxKillTU)>0                 % if killing happens, then update  
        L(TUcells(auxKillTU)) = false;  % Remove the killed tumor cell from the grid
        Lt(TUcells(auxKillTU)) = false;  % Additionally remove the killed tumor cell from the tumor cell grid
        TUcells(auxKillTU) = [];            % remove from stack of tumor cells
        %Remove properties of the killed tumor cells
        TUprop.isStem(auxKillTU) = []; % remove stemness property
        TUprop.Pcap(auxKillTU) = [];  % remove Pmax property
        TUprop.mut(auxKillTU) = []; %Remove the mutation property
        TUprop.EMT(auxKillTU) = []; %Remove the EMT property 
        Mprop.Kcap(auxKillM) = Mprop.Kcap(auxKillM)-1; % Decrease the killing capacity of macrophages 
        Mprop.engaged(auxKillM) = engagementDuration; % killers become engaged (it takes some time before they can kill again)
    else
        auxKillTU = 0;
    end
    else
        auxKillTU = 0;
    end % Reset killing filter
else
    auxKillTU = 0;
end % Reset candidate filter


end