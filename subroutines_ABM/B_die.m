function [L, Bcells] = B_die(L,Bcells,Bpdeath,nh)
%Function in which basement membrane performs its actions
%As described in the paper, there is currently no chance of the basement
%membrane dying spontaneously (Bpdeath=0), but if you were to prescribe
%actions to the basement membrane, you can do it using this function. 

m = getAdjacent(L,Bcells,nh);
%Determine neighbourhood of the basement membrane agents 
%If an agents wants to migrate or proliferate for example, it requires an
%empty spot in its neighborhood. 

D = m.randI2 <= (Bpdeath); %Generate random numbers to compare to the pdeath thresold
%D contains all agents that will perform the death action

if ~isempty(D) % only go through loop if there are cells that will perform this action 
    L(Bcells(m.S2(D))) = false;    % remove basement membrane agent from grid
    Bcells(m.S2(D)) = [];          % remove from stack of basement membrane agents
end