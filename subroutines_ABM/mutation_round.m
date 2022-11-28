function [TUcells,TUprop] = mutation_round(TUcells,TUprop,TUpmut,TUantig,L,nh,TUmutmax)
%Function that describes the mutating of the luminal cells. It also adds the proliferative advantage obtained with mutations and returns an updated
%list of luminal cells. 

%CREATE MASKS FOR ADJACENT POSITIONS
m = getAdjacent(L,TUcells,nh);

%Include added chance of gaining a mutation for cells that already mutated
factor = 0.05*3;
MapMut = TUprop.mut*factor + TUpmut;

mut = m.randI2 <= MapMut & (TUprop.mut(m.S2)<TUmutmax); %Select the cells that will mutate this round 
factor=50*3;
if sum(mut)>0
    TUprop.mut(m.S2(mut)) = TUprop.mut(m.S2(mut))+TUantig; %Update the mutation property of mutated cells 
    TUprop.Pcap(m.S2(mut))= TUprop.Pcap(m.S2(mut))+TUantig*factor; %Update the altered proliferation capacity of the newly mutated cells 
end

    