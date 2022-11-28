function [Bcells, L,Lb] = TU_kill(Bcells, L,K,m,ChtaxMap,Lb,TUcells)
%Function that allows tumor cells that have gained enough mutations to kill basement membrane agents (to allow for invasion). Function returns an
%updated list of basement membrane agents and an updated grid. 

K2 = find(K==1); %Select the tumor cells that will kill a basement membrane agent
N=125;

for l=1:numel(K2) %only for the cells that will actually do something
    curID = K2(l); % number within stack of currently acting cell
    nh = [TUcells(curID)-1,TUcells(curID)+1,TUcells(curID)+N,TUcells(curID)-N]; %To only select 4 neighbors (not Moore's neighborhood)
    nh2 = nh(nh>0); % erasing places that were previously occupied by a tumor cell (cells cannot move to occupied spaces)
    chem = ChtaxMap(nh2(:));
    nh3 = nh2(chem==0); %select the spots in the neighborhood that are occupied by a basement membrane agent
    if ~isempty(nh3)
        L(nh3(1)) = false; %Remove the basement membrane agent from the grid
        Lb(nh3(1)) = false; %Remove the basement membrane agent from the BM grid 
        id = find(Bcells==nh3(1)); 
        Bcells(id) = []; %Remove the basement membrane agent from the list of basement membrane agents
    end
end

end
