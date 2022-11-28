function [mySystem,TUcellNo,FcellNo,finalImage,TSC,Tm,Tm2,MutCells,BM] = growTumor(mySystem,cnst)
%Performs the agent based modeling.
%Inputs are two structures; the parameters and some global constants defined in getSystemParams.
%Ouptut is the system at final iteration, the number of tumor cells over time (other selected parameters at the bottom of this script) and the final video. 

%BUILD ENVIRONMENT
%Throw parameters to workspace:
cellfun(@(f) evalin('caller',[f ' = mySystem.params.' f ';']), fieldnames(mySystem.params));
cellfun(@(f) evalin('caller',[f ' = mySystem.grid.' f ';']), fieldnames(mySystem.grid));

%INITIALIZE SYSTEM
if cnst.newSystem
    %Create empty system:
    L=false(N,M);                                                      %initialize grid
    L(1,:) = true; L(:,1) = true; L(end,:) = true; L(:,end) = true;     %set grid boundary to occupied
    L2=false(N,M);                                                      %initialize grid
    L2(1,:) = true; L2(:,1) = true; L2(end,:) = true; L2(:,end) = true;     %set grid boundary to occupied

    %Create lists for placements of the different agents
    placesF = []; %Fibroblasts
    placesC = []; %Basal cells
    placesB = []; %Basement membrane
    placesTU = []; %luminal cells
    duct = []; %Lumen inside the acinus (first half)
    duct2 = []; %Lumen inside the acinus (second half)
    degd=[]; %Possible empty spaces in the stroma resulting from fibroblasts
    degd2 = []; %Possible empty spaces in the stroma resulting from tumor cells
    
    %Create prostatic acinus 
    Fcellcenter2=N*round(M/2)-round(N/2); %Define center of the acinus 
    for q = 1:mySystem.params.Number %Create the number of acini as defined in the parameter file 
        Fcellcenter=Fcellcenter2(q);
        %Placement of the cell is done in four quadrants, defined below 
        centr1 = Fcellcenter + Fradius;
        centr2 = Fcellcenter - Fradius;
        centr3 = Fcellcenter - (Fradius)*N;
        centr4 = Fcellcenter + (Fradius)*N;
        %Build up the acinus working from the most outward layer of the
        %acinus to the most inward one : placing the center agents
        placesF = [placesF,centr1,centr2,centr3,centr4]; %Outer layer of fibroblasts surrounding the acinus
        placesB = [placesB,centr1-1,centr2+1,centr3+N,centr4-N]; %Layer of basement membrane
        placesC = [placesC,centr1-2,centr2+2,centr3+2*N,centr4-2*N]; %layer of basal cells 
        placesTU = [placesTU,centr1-3,centr1-4,centr1-5,centr2+3,centr2+4,centr2+5,centr3+3*N,centr3+4*N,centr3+5*N,centr4-3*N,centr4-4*N,centr4-5*N]; %layer of tumor cells (2 cell layers thick)
        placesF = [placesF,centr1+1,centr2-1,centr3-N,centr4+N]; %Second layer of fibroblasts at the outside of the acinus
        duct = [duct,centr2+6:centr1-6]; %empty spaces inside of the acinus 

        for i = 1:Fradius
            number=round(sqrt((Fradius^2)-(i^2))); %Determining delta y for each of the quadrants
            number = Fradius - number ;
            %Start moving from the center (starting point of the quadrant
            %to the next quadrant) - for each quadrant this requires a
            %different movement, described below:
            coord1= centr1 + (i*N) - number; 
            coord2 = centr1 - (i*N) - number;
            coord3 = centr2 + (i*N) + number;
            coord4 = centr2 - (i*N) + number;
            placesF=[placesF, coord1,coord2,coord3,coord4]; %Add fibroblasts in the outer layers
            placesF=[placesF, coord1+1,coord2+1,coord3-1,coord4-1]; %Add fibroblasts in the outer layers
            if i < Fradius -1 %Move inwards towards the basement membrane layer
                placesB=[placesB,coord1-1,coord2-1,coord3+1,coord4+1]; 
            end
            if i < Fradius-2 %Move inwards towards the basal cell layer
                placesC=[placesC,coord1-2,coord2-2,coord3+2,coord4+2];
            end
            if i < Fradius -3 %Move inwards towards the luminal cell layer (1)
                placesTU = [placesTU,coord1-3,coord2-3,coord3+3,coord4+3];
            end
            if i < Fradius -4 %Move inwards towards the second luminal cell layer 
                placesTU = [placesTU, coord1-4,coord2-4,coord3+4,coord4+4];
            end
            if i < Fradius -5
                placesTU =[placesTU, coord1-5,coord2-5,coord3+5,coord4+5];
            end
            if i < Fradius -5
                Fr = (i*2); %Determine remaining space 
                duct = [duct,coord1-6:coord3+6,coord4+6:coord2-6]; %Set it to empty spaces aka duct
                duct2 = [duct2,coord1-6+(Fr*N):coord3+6+(Fr*N),coord4+6+(Fr*N):coord2-6+(Fr*N)]; 
            end
        end

        for j = 1:Fradius
            number=round(sqrt((Fradius^2)-(j^2)));
            number = Fradius - number;
            %Do the other half of the movements from the center (starting point of the quadrant
            %to the next quadrant) - for each quadrant this requires a
            %different movement, described below:
            coord1 = centr3 - j +(N*number);
            coord2 = centr3 + j +(N*number);
            coord3 = centr4 - j -(N*number);
            coord4 = centr4 + j - (N*number);
            placesF=[placesF,coord1,coord2,coord3,coord4]; %Add fibroblasts to the outer two layers
            placesF=[placesF,coord1-N,coord2-N,coord3+N,coord4+N];
            if j< Fradius -1 %Move inwards towards the basement membrane
                placesB=[placesB,coord1+N,coord2+N,coord3-N,coord4-N];
            end
            if j< Fradius - 2 %Move inwards towards the basal cell layer
                placesC=[placesC, coord1+2*N,coord2+2*N,coord3-2*N,coord4-2*N];
            end
            if j< Fradius - 3 %Move inwards towards the luminal cell layers
                placesTU=[placesTU,coord1+3*N,coord2+3*N,coord3-3*N,coord4-3*N];
            end
            if j< Fradius -4
                placesTU=[placesTU,coord1+4*N,coord2+4*N,coord3-4*N,coord4-4*N];
            end
            if j< Fradius -5
                placesTU=[placesTU,coord1+5*N,coord2+5*N,coord3-5*N,coord4-5*N];
            end
        end

        %Make sure there are no double numbers within each vector or
        %between vectors (there can only be one agent in each spot)
        placesF=unique(placesF);

        [x,y1]=intersect(placesB,placesF); %Compare fibroblast layer and basement membrane layer
        placesB(y1)=[]; %Remove possible double spots from the basement membrane layer

        [x,y2]=intersect(placesC,placesB); %Compare basement membrane with basal cell layer
        placesC(y2)=[]; %Remove possible double spots from the basal cell layer

        [x,y3] = intersect(placesTU,placesC); %Compare basal cell layer with tumor cell layer
        placesTU(y3)=[]; %Remove possible double spots from the tumor cell layer 

        [xc,yc]=intersect(placesTU,placesB); %Compare basement membrane with tumor cell layer
        placesTU(yc)=[]; %Remove possible double spots from the tumor cell layer

        placesB2=unique(placesB);
        placesC2=unique(placesC);
        placesTU2=unique(placesTU);

        coord = Fcellcenter+Fradius+10;
        %If the desired shape is an ellips and not a circle, adapt the
        %shape to an ellips (two different 'radii')
        if mySystem.params.Shape == 2
            placesB2(placesB2<coord)=placesB2(placesB2<coord)-N*Fellips;
            placesB2(placesB2>coord)=placesB2(placesB2>coord)+N*Fellips;
            placesC2(placesC2<coord)=placesC2(placesC2<coord)-N*Fellips;
            placesC2(placesC2>coord)=placesC2(placesC2>coord)+N*Fellips;
            placesTU2(placesTU2<coord)=placesTU2(placesTU2<coord)-N*Fellips;
            placesTU2(placesTU2>coord)=placesTU2(placesTU2>coord)+N*Fellips;
            placesF(placesF<coord)=placesF(placesF<coord)-N*Fellips;
            placesF(placesF>coord)=placesF(placesF>coord)+N*Fellips;
            duct(duct<coord)=duct(duct<coord)-N*Fellips;
            duct2(duct2>coord)=duct2(duct2>coord)+N*Fellips;
            placesF=[placesF,((Fcellcenter-Fradius)-N*(Fellips-1):N:(Fcellcenter-Fradius)+N*Fellips)];
            placesF=[placesF,((Fcellcenter+Fradius)-N*(Fellips-1):N:(Fcellcenter+Fradius)+N*Fellips)];
            placesF=[placesF,((Fcellcenter-Fradius-1)-N*(Fellips-1):N:(Fcellcenter-Fradius-1)+N*Fellips)];
            placesF=[placesF,((Fcellcenter+Fradius+1)-N*(Fellips-1):N:(Fcellcenter+Fradius+1)+N*Fellips)];
            placesB2=[placesB2,((Fcellcenter-Fradius+1)-N*(Fellips-1):N:(Fcellcenter-Fradius+1)+N*Fellips)];
            placesB2=[placesB2,((Fcellcenter+Fradius-1)-N*(Fellips-1):N:(Fcellcenter+Fradius-1)+N*Fellips)];
            placesC2=[placesC2,((Fcellcenter-Fradius+2)-N*(Fellips-1):N:(Fcellcenter-Fradius+2)+N*Fellips)];
            placesC2=[placesC2,((Fcellcenter+Fradius-2)-N*(Fellips-1):N:(Fcellcenter+Fradius-2)+N*Fellips)];
            placesTU2=[placesTU2,((Fcellcenter-Fradius+3)-N*(Fellips-1):N:(Fcellcenter-Fradius+3)+N*Fellips)];
            placesTU2=[placesTU2,((Fcellcenter+Fradius-3)-N*(Fellips-1):N:(Fcellcenter+Fradius-3)+N*Fellips)];
            placesTU2=[placesTU2,((Fcellcenter-Fradius+4)-N*(Fellips-1):N:(Fcellcenter-Fradius+4)+N*Fellips)];
            placesTU2=[placesTU2,((Fcellcenter+Fradius-4)-N*(Fellips-1):N:(Fcellcenter+Fradius-4)+N*Fellips)];
            placesTU2=[placesTU2,((Fcellcenter-Fradius+5)-N*(Fellips-1):N:(Fcellcenter-Fradius+5)+N*Fellips)];
            placesTU2=[placesTU2,((Fcellcenter+Fradius-5)-N*(Fellips-1):N:(Fcellcenter+Fradius-5)+N*Fellips)];
            for k = 1:Fellips
                duct=[duct,(Fcellcenter-Fradius+6)-N*(Fellips-k):(Fcellcenter+Fradius-6)-N*(Fellips-k)];
                duct2=[duct2,(Fcellcenter-Fradius+6)+N*(Fellips+1-k):(Fcellcenter+Fradius-6)+N*(Fellips+1-k)];
            end
        end
    end
    
    placesTU3=placesTU2<N; %Make sure no of the tumor cells fall out of the grid 
    placesTU2(placesTU3)=[]; 
    
    %Create luminal cell agent
    TUcells = int32(placesTU2); %Add all the placements of the luminal cells to the list
    TUcellNo = length(TUcells); %Determine the number of luminal cells 
    TUprop.isStem = false(1,TUcellNo);  %Set all luminal cells to not being a stem cell initially
    randng = (0.5/TUcellNo)+rand(1,round(0.3*TUcellNo))*(1-(0.5/TUcellNo)); %Set the defined number of luminal cells to being a stem cell
    for i=1:round(03*TUcellNo) %Choose placement of stem cells randomly
        TUprop.isStem(round(randng*TUcellNo))=true;
    end
    %Set luminal cell properties
    TUprop.Pcap = uint8(TUpmax*ones(1,TUcellNo)); %Maximum proliferation capacity for each cell
    TUprop.mut = zeros(1,TUcellNo); %Amount of mutations for each cell (initially all 0)
    TUprop.EMT = zeros(1,TUcellNo); %Whether or not a cell has gone through EMT (initially none have done this)

    placesF3=placesF<N; %Make sure no of the fibroblasts fall out of the grid 
    placesF(placesF3)=[]; 
    
    %Create fibroblast agent
    Fcells=int32(placesF); %Add all the placements of the fibroblasts to the list
    FcellNo=length(Fcells);%Determine number of fibroblasts
    Fprop.Pcap = uint8(repmat(Fpmax,1,FcellNo)); %Set proliferation capacity of fibroblasts as a property 
    For = placesF;
    
    C = (find(L==false))'; %Find empty spots on the grid 
    CFcells= int32(C(randperm(numel(C),CFcellNo)));    %randomly allocate positions for CAFs (initial number is generally 0, so usually an empty vector)
    CFprop.Pcap = uint8(repmat(CFpmax,1,CFcellNo));    %add properties: max killing capacity
    
    placesB3 = placesB2<N; %Make sure no of the basement membrane falls out of the grid
    placesB2(placesB3)=[]; 
    Bcells = int32(placesB2); %Add all the placements of the basement membrane to the list
    BcellNo=length(Bcells); %Determine the number of basement membrane agents
    
    placesC3=placesC2<N; %Make sure no of the basal cells fall out of the grid 
    placesC2(placesC3)=[];
    Ccells =int32(placesC2); %Add all the placements of the basal cells to the list
    CcellNo=length(Ccells); %Determine the number of basal cell agents
    Cprop.Pcap = uint8(repmat(Cpmax,1,CcellNo)); %Set the property of the basal cells (maximum proliferation capacity
    
    %Randomly add defined number of macrophages type 1 and 2 to grid
    C = (find(L==false))' ;                            %check for free spots on the grid
    M1cells= int32(C(randperm(numel(C),M1cellNo)));    %randomly allocate positions for initially present M1 macrophages
    M1prop.Kcap = uint8(repmat(M1kmax,1,M1cellNo));    %add properties: max killing capacity
    M1prop.engaged = uint8(zeros(1,M1cellNo));         %add properties: engagement in killing

    C = (find(L==false))';                             %check for free spots on the grid
    M2cells= int32(C(randperm(numel(C),M2cellNo)));    %randomly allocate positions for M1s
    M2prop.Kcap = uint8(repmat(M2kmax,1,M2cellNo));    %add properties: max killing capacity
    M2prop.engaged = uint8(zeros(1,M2cellNo));         %add properties: engagement in killing
    
    %Randomly place fibroblasts in the stroma
    Fcellsst = int32(C(randperm(numel(C),round(0.01*(M*N)))));
    %Make sure none of the fibroblasts overlap with the places occupied by
    %the acinus
    [xa,ya]=intersect(Fcellsst,placesC2); Fcellsst(ya)=[];
    [xb,yb]=intersect(Fcellsst,placesTU2); Fcellsst(yb)=[];
    [xc,yc]=intersect(Fcellsst,placesB2); Fcellsst(yc)=[];
    [xd,yd]=intersect(Fcellsst,duct); Fcellsst(yd)=[];
    [xe,ye]=intersect(Fcellsst,duct2); Fcellsst(ye)=[]; 
        
    Fcells=[Fcells, Fcellsst]; %Add the fibroblasts in the stroma to the list of fibroblasts
    FcellNo=length(Fcells); %Update the number of fibroblasts 
    Fprop.Pcap = [Fprop.Pcap, uint8(repmat(Fpmax,1,length(Fcellsst)))]; %Update the fibroblast properties (add for the newly placed agents)
    
    if cnst.video
        figure()
    end
        
else
    %Use existing system to grow furhter
    cellfun(@(f) evalin('caller',[f ' = mySystem.TU.' f ';']), fieldnames(mySystem.TU));
    cellfun(@(f) evalin('caller',[f ' = mySystem.F.' f ';']), fieldnames(mySystem.F));
    cellfun(@(f) evalin('caller',[f ' = mySystem.B.' f ';']), fieldnames(mySystem.B));
    cellfun(@(f) evalin('caller',[f ' = mySystem.C.' f ';']), fieldnames(mySystem.C));
    cellfun(@(f) evalin('caller',[f ' = mySystem.M1.' f ';']), fieldnames(mySystem.M1)); 
    cellfun(@(f) evalin('caller',[f ' = mySystem.M2.' f ';']), fieldnames(mySystem.M2));
    cellfun(@(f) evalin('caller',[f ' = mySystem.CF.' f ';']), fieldnames(mySystem.CF));
    cellfun(@(f) evalin('caller',[f ' = mySystem.grid.' f ';']), fieldnames(mySystem.grid));
end

%Place all the cells on the grid    
L(TUcells)=true;                            %place tumor cells on grid
Lt = false(size(L));                        %reset tumor grid
Lt(TUcells) = true;                         %update tumor grid

L([TUcells,Fcells])=true;                   %add fibroblasts to grid
Lf = false(size(L));                        %reset fibroblast grid
Lf(Fcells) = true;                          %update fibroblast grid

L([TUcells,Fcells,Bcells])=true;            %Add basement membrane to the grid
Lb = false(size(L));                        %Reset basement membrane grid
Lb(Bcells) = true;                          %Update basement membrane grid

L([TUcells,Fcells,Bcells,Ccells])=true;     %Place basal cells on the grid
Lc = false(size(L));                        %Reset basal cell grid
Lc(Ccells) = true;                          %Update basal cell grid

L([TUcells, Fcells, Bcells, Ccells, M1cells])=true;         %add M1 cells to grid
Lm1 = false(size(L));                       %Reset M1 grid
Lm1(M1cells) = true;                        %Update M1 grid

L([TUcells, Fcells, Bcells, Ccells, M1cells, M2cells])=true;%add M2 cells to grid
Lm2 = false(size(L));                       %Reset M2 grid
Lm2(M2cells) = true;                        %Update M2 grid

L([TUcells, Fcells, Bcells, Ccells, M1cells, M2cells, CFcells]) = true; %Add CAFs to the grid
Lcf = false(size(L));                       %Reset CAF grid
Lcf(CFcells) = true;                        %Update CAF grid

%seeding of stroma 
Ls = true(size(L));                         %Create full grid
Ls([TUcells, Fcells, Bcells, Ccells, M1cells, M2cells, CFcells, duct, duct2,253]) = false; %Remove places occupied by the other cells
original = [TUcells,Fcells,Bcells,Ccells,M1cells,M2cells,CFcells,duct,duct2];       %Save starting configuration
all = 1:(125*125);
Lstin = ismember(all,original);
all(Lstin)=[]; %All spots that are already occupied are removed
L([TUcells, Fcells, Bcells, Ccells, M1cells, M2cells, CFcells, all]) = true; %Spots that were not occcupied (outside the acinus) are now occupied by stroma

%seeding of stroma - immune
Ls2 = true(size(L)); %Create a second variable for the stroma to track permeability for immune cells
Ls2([TUcells, Fcells, Bcells, Ccells, M1cells, M2cells, CFcells, duct, duct2,253]) = false; %Only occupy spaces that were not already taken by agents 
    
%Set auxiliary variables:
nh.aux = int32([-N-1 -N -N+1 -1 1 N-1 N N+1])'; % indices the Moore neighbourhood of an agent
nh.Pms = perms(uint8(1:8))';                % permutations of adjacent positions
nh.nP = size(nh.Pms,2);                     % number of possible permutations
fcount = 0;                                 %frame counter for video
finalImage = [];                            %set empty resulting image

%second neighborhood (two times the Moore neighbourhood)
nh2 = int32([-2*N-2 -2*N-1 -2*N -2*N+1 -2*N+2 -N-2 -2 N-2 -N+2 2 N+2 2*N-2 2*N-1 2*N 2*N+1 2*N+2 -N-1 -N -N+1 -1 1 N-1 N N+1])'; 

%ITERATIONS THROUGH NUMBER OF TIME STEPS
%Create empty lists for variables that you want to save
Tm2 = 0; 
killstt=[0,0];
BM = [];
for i = 1:cnst.nSteps                     %iterate through time steps
    
%     if (sum(M2cells)~=0 && M2TUadd > 0) %Check wether there are M2Cells present and if they are tumor promoting
            %Whether or not you want to use this method depends on from how
            %far M2 macrophages affect tumor cell proliferation
%         TUpprol = mySystem.params.TUpprol+M2TUadd; %
%     else %if not tumor promoting M2cells present
%         TUpprol = mySystem.params.TUpprol;
%     end
    if killstt(i+1)>killstt(i) %After having found and killed the first mutated cell, M1 macrophage influx is increased (responding to the 'bad cell')
        M1influxProb=M1influxProb*3;
    end

    if (Tm2(i)/length(TUcells))>0.3 %Once more than 30% of the luminal cells has mutated, basal cell equilibrium is disturbed (they can now die and their spot can be taken up by mutated cells)
        Cpdeath=0.01;
    end
    

    %Create chemotaxis maps
    ChtaxMapF = double(bwdist(Lf,'euclidean')); %Fibroblasts 
    ChtaxMapC = double(bwdist(Lc,'euclidean')); %Basal cells 
    ChtaxMapB = double(bwdist(Lb,'euclidean')); %Basement membrane
    ChtaxMapM1 = double(bwdist(Lm1,'euclidean')); %M1 macrophages
    ChtaxMapM2 = double(bwdist(Lm2,'euclidean')); %M2 macrophages 
    ChtaxMapCF = double(bwdist(Lcf, 'euclidean')); %Fibroblasts 
    
    degd3 = unique([degd,degd2]); %Find spaces in the stroma that were affected by tumor cells/fibroblasts and possibly made space for other cells (permeability)
    degd3(degd3<N)=[]; %Remove all spaces outside of the grid 
    degd3(degd3>N*M)=[];
    Ls = true(size(L));
    Ls(unique([original,TUcells, Fcells, Bcells, Ccells, M1cells, M2cells, CFcells, duct, duct2,253])) = false; %Make sure to not remove cells from the acinus (only stroma)
    rndnrs = rand(size(L)); %obtain random numbers
    mapp = ones(size(rndnrs));
    mapp(degd3)=100*2;
    newrnd=rndnrs./mapp;%Make the random numbers of the selected spots (degradation by other cells) extra small, so that these will definitely be free
 
    rm = newrnd<stromaPerm; %Compare random numbers to probability of stromal permability 
    taken = sort(unique([original,TUcells,Fcells,Bcells,Ccells,M1cells,M2cells,CFcells,duct,duct2]));
    taken2 = zeros(size(L));
    taken2(taken)=1; %Set all values that were originally taken (acinus) in the starting geometry to 1

    taken2(rm==1)=1; %additionally, set all spaces in the stroma that yielded probabilities below the stromaPerm to 1
    all = 1:(125*125);
    all(taken2==1) =[]; %No stroma will be placed at the original geometry (acinus) and at the recently freed up spaces

    Ls(newrnd<stromaPerm) = false; %Remove the freed up spaces from the stromal vector 
    ChtaxMapS = double(bwdist(Ls, 'euclidean')); %Update the stromal chemotaxis map
    L(Ls) = true; 
    
    %Luminal cell round
    L([TUcells,M1cells,M2cells,Fcells,Bcells, Ccells]) = true;      %ensure that all cells are on the grid
    [TUcells,TUprop] = shuffleCells(TUcells,TUprop); %shuffle tumor cells
    [L, TUcells, TUprop, TUdeath,Bcells,Lb] = TU_go_grow_die(L, nh, TUcells, TUprop, TUpprol, TUpmig, TUpdeath, TUps, TUrwalk, ChtaxMapF,ChtaxMapC,TUwidth,ChtaxMapB,TUpkill,Bcells,Lb,TUthrshBM,For,BM); %Update grid with growing tumor
    Lt = false(size(L));                    %reset luminal cell grid
    Lt(TUcells) = true;                     %update luminal cell grid
    %Create luminal cell chmotaxis map
    ChtaxMapTU = double(bwdist(Lt,'euclidean'));
    
    %Interaction round of the luminal cells with the other cells
    %(additional actions described in the paper)
    [L, M1cells, M1prop, M2cells, M2prop,degd] = TU_interaction_round(L, nh, TUcells, TUpMdiff, TUthrshM, M1cells, M1prop, M2cells, M2prop, ChtaxMapM1,TUprop,TUpMMP,TUthrshMMP,degd,ChtaxMapS);
     
    %Fibroblast round
    [Fcells,Fprop] = shuffleCells(Fcells,Fprop);%shuffle fibroblasts
    [L, Fcells, Fprop,CFcells,CFprop] = F_go_grow_die(L, nh, Fcells, Fprop, Fpprol, Fpmig, Fpdeath, Frwalk, ChtaxMapTU,CFcells,CFprop,Fdiff);
    
    %CAF round
    [CFcells, CFprop] = shuffleCells(CFcells, CFprop); %Shuffle CAFs
    [L, CFcells, CFprop,TUcells, TUprop,degd2,M1cells, M1prop, M2cells, M2prop] = CAF_go_grow_die(L, nh, CFcells, CFprop, CFpprol, CFpmig, CFpdeath, CFrwalk, ChtaxMapTU,TUcells,TUprop,CFemt,ChtaxMapCF,CFprom,ChtaxMapS,CFmmp,ChtaxMapM1,CFmdiff,M1cells, M1prop, M2cells, M2prop,BM);
    Lf = false(size(L));  %Reset fibroblast grid
    Lf(Fcells) = true;  
    Lcf = false(size(L)); %Reset CAF grid
    Lcf(CFcells) = true;
    
%     %Basement membrane round (now a quiescent agent)
%You should create B.prop if you want to use ShuffleCells
%     [Bcells] = shuffleCells(Bcells); %Shuffle basement membrane
%     [L, Bcells] = B_die(L,Bcells,Bpdeath,nh); 
%     Lb = false(size(L)); %Update basement membrane grid
%     Lb(Bcells) = true;
%     ChtaxMapB = double(bwdist(Lb,'euclidean'));
    
    %Basal cell round
    [Ccells, Cprop] = shuffleCells(Ccells,Cprop); %Shuffle basal cells
    [L, Ccells, Cprop] = C_go_grow_die(L, nh, Ccells, Cprop, Cpprol, Cpmig, Cpdeath, Crwalk, ChtaxMapF, ChtaxMapB);
    Lc = false(size(L)); %Update basal cell grid 
    Lc(Ccells) = true;
    
    %Update grid for immune cell permeability of the stroma
    Ls2 = true(size(L2));
    Ls2(unique([original,TUcells, Fcells, Bcells, Ccells, M1cells, M2cells, CFcells, duct, duct2,253])) = false;
    rndnrs = rand(sum(Ls2(:)),1);
    L2(Ls2) = rndnrs>stromaImPerm; %All stroma that has a value below the StromaImPerm is free for immune cells to move to 
    L2(1,:) = true; L2(:,1) = true; L2(end,:) = true; L2(:,end) = true;
    
    %M1 and M2 round
    [L,M1cells,M1prop,TUcells,TUprop,Lt,killst] = mCellRound(L2,Lt,M1cells,TUcells,M1prop,TUprop,M1influxProb,M1influxRate,M1kmax,M1speed,M1pmig,M1rwalk,M1pkill,M1pdeath,M1engagementDuration,ChtaxMapTU,nh,M2cells,Fcells,nh2);
    L([TUcells,M1cells,M2cells,Fcells,Ccells,Bcells]) = true;      %ensure that all cells are on the grid
    Lm1 = false(size(L)); %Update M1 cell grid 
    Lm1(M1cells) = true; 
    killstt=[killstt,killst]; %Keep track of killing of mutated cells (to possibly enhance M1 cell influx)

    [L,M2cells,M2prop,TUcells,TUprop,Lt] = mCellRound2(L2,Lt,M2cells,TUcells,M2prop,TUprop,M2influxProb,M2influxRate,M2kmax,M2speed,M2pmig,M2rwalk,M2pkill,M2pdeath,M2engagementDuration,ChtaxMapTU,nh2,M1cells,Fcells,nh,M2TUadd,ChtaxMapM2,M2emt,BM);
    L([TUcells,M1cells,M2cells,Fcells,Ccells,Bcells]) = true;      %ensure that all cells are on the grid
    Lm2 = false(size(L)); %Update M2 cell grid
    Lm2(M2cells) = true;
    
    %Mutation round
    L([TUcells,M1cells,M2cells,Fcells,Ccells,Bcells,CFcells]) = true;      %ensure that all cells are on the grid
    [TUcells,TUprop] = mutation_round(TUcells,TUprop,TUpmut,TUantig,L,nh,TUmutmax);
    
    %Update the system 
    mySystem=updateSystem(mySystem,TUcells,Fcells,Bcells,Ccells,M1cells,M2cells,CFcells,TUprop,M1prop,M2prop,Fprop,Cprop,CFprop,L,Lt,Ls,i);

    %Store variables each iteration 
    TUcellNo(i+1,1)= length(TUcells); %Number of luminal cells
    FcellNo(i+1,1) = length(Fcells); %Number of fibroblasts
    TSC(i+1,1)=length(find(TUprop.isStem==1));
    %TSCmut(i+1,1) = length(find(TUprop.isStem==1 & TUprop.mut>0));
    %%Percentage of mutated stem cells
    Tm(i+1,1)=mean(TUprop.mut);
    Tm2(i+1,1)=length(find(TUprop.mut>0)); %number of mutated cells
    BM(i+1,1)=length(Bcells);   %number of basement membrane agents
    MutCells(i+1,1) = sum(TUprop.mut)/length(find(TUprop.mut>0)); %Average amount of mutations per mutated cell
   
%DRAW IF OPTION IS SELECTED
fcount=fcount+1;
if cnst.video
    visualizeSystem(TUprop,mySystem,all);
    drawnow, currFrame = getframe(gca);
    finalImage{fcount} = currFrame.cdata;
else
    finalImage{fcount} = false;    
end
end 
    

end

