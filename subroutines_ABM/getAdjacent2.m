function m = getAdjacent2(L,MYcells,nh2)
%Creates masks for the inputted list of cells (e.g. neighborhood space) and
%creates lists of random numbers for all cells in the list and a separate one for the agents that have a free spot in their neighborhood. The created
%masks are provided as output of this function (m.). Similar function to getAdjacent, but this time for a two times Moore neighborhood

for k = 1:length(MYcells)
    vect=nh2'+MYcells(k);
    Delete=vect<1 | vect>15625; %Select spaces that should not be included in the neighborhoods (at the boundaries or outside of the grid)
    vect(Delete)=[]; %Remove spaces that should not be included 
    if length(vect)<24
        vect = [vect,double(vect(length(vect)))*ones(1,24-length(vect))]; % Since the Moore neighborhood is twice as large, it can be that some cells have fewer neighbors (fill the list with 'dummies' so that all vectors are of same length)
    end
    m.S(:,k)=vect'; %save the neighborhood for the selected cell 
end
m.S2 = 1:length(MYcells);   % List of all agents 
m.S(L(m.S)) = 0; 			% setting occupied grid cells to false
m.indxF = find(any(m.S)); 	% selecting agents with at least one free spot
m.nC = length(m.indxF); 	% number of agents with free spot
m.randI = rand(1,m.nC); 	% initialize random number vector for agents with at least one free spot in their neighborhood
m.nC2 = length(MYcells);    % Determining amount of agents
m.randI2 = rand(1,m.nC2);   % Initialize random number vector for all agents
end