function m = getAdjacent(L,MYcells,nh)
%Creates masks for the inputted list of cells (e.g. neighborhood space) and
%creates lists of random numbers for all cells in the list and a separate one for the agents that have a free spot in their neighborhood. The created
%masks are provided as output of this function (m.)

m.S = bsxfun(@plus,MYcells,nh.aux(nh.Pms(:,randi(nh.nP,1,length(MYcells)))));
m.S3 = m.S;                 % Copying m.S (to save an unaltered version of it)
m.S2 = 1:length(MYcells);   % List of agents that have no free spot in their neighborhood (e.g. for death an empty space is not required)
m.S(L(m.S)) = 0; 			% setting occupied grid cells to false
m.indxF = find(any(m.S)); 	% selecting agents with at least one free spot
m.nC = length(m.indxF); 	% number of agents with free spot
m.randI = rand(1,m.nC); 	% initialize random number vector for all cells with an empty space in their neighborhood
m.nC2 = length(MYcells);    % Determine total amount of agents
m.randI2 = rand(1,m.nC2);   % Create a list of random numbers for all agents (also the ones with no free spot in their neighborhood)
end