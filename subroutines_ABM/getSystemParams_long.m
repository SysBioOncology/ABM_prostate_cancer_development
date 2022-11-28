function [mySystem, cnst]=getSystemParams_long(space)
%Sets the parameters needed for modeling the ABM (using timesteps of 12 hours)

mySystem.params.initialSeed = 1;   % initial random seed, default 1
mySystem.params.stromaPerm = 0.001; % stroma permeability in general (chance that an agent can take up a space that was stroma initially)
mySystem.params.stromaImPerm = 0.8; %stroma permeability for immune cells 
mySystem.params.Shape = 1; %Setting the shape of the acinus: 2 indicating ellips, 1 indicating circle 
mySystem.params.Number = 1; %Indicating the number of acini that you want 
mySystem.params.Fradius = 0; %Radius of the acinus
mySystem.params.Fellips= 0;  %Second 'radius' of the acinus, if you are plotting an ellips

%DEFINE TUMOR CELL PARAMETERS (luminal cells)
mySystem.params.TUpprol = 0.3054;  %probability of proliferation of luminal cells
mySystem.params.TUpmig = 0.1167*3;     %probability of migration
mySystem.params.TUpdeath = 0.00284;  %probability of death
mySystem.params.TUrwalk = 0.5;    %random influence on movement
mySystem.params.TUpmax = 4;   %initial proliferation capacity
mySystem.params.TUps = 0.2967; %probability of symmetric division (for stem cells)
mySystem.params.TUcellNo = 1500;      %Number of initial luminal cells, only used when randomly seeding luminal cells (otherwise it plots the shape of the acinus)
mySystem.params.TUpmut = 0.00015; %probability of acquiring a mutation  
mySystem.params.TUantig = 0.005; %Added 'antigenicity' per mutation 
mySystem.params.TUmutmax = 1; %maximum antigenicity a cell can have (i.e. TUmutmax/Tuantig = amount of mutations one cell can acquire)
mySystem.params.TUwidth = 2; %Width of the layer of initial healthy luminal cells
mySystem.params.TUpkill = 0.3; %Chance of breaking down Basement membrane upon reaching the amount of required mutations
mySystem.params.TUthrshBM = mySystem.params.TUantig*6; %Threshold of accumulated mutations value for breaking down BM
mySystem.params.TUpMdiff = 0.8; %Chance that a tumor cell affects macrophage differentiation (via cytokines) from M1 to M2
mySystem.params.TUthrshM = mySystem.params.TUantig*14; %Needed threshold of accumulated mutations for differentiation effect
mySystem.params.TUpMMP = 0; %Chance that tumor cells create a path in the stroma to migrate to/through
mySystem.params.TUthrshMMP = mySystem.params.TUantig*5; %Required amount of mutations before being able to break down stroma

%DEFINE MACROPHAGE TYPE 1 PARAMETERS
mySystem.params.M1kmax = 11;        %killing capacity 
mySystem.params.M1pkill = 0.0918; %probability of killing 
mySystem.params.M1pmig = 0.8001; %probability of migration
mySystem.params.M1pdeath = 0.0147; %probability of death
mySystem.params.M1rwalk = 0.8;      %random influence on movement
mySystem.params.M1speed = 40;       %speed of movement (amount of times a cell can move per timestep)
mySystem.params.M1influxProb = 0.08; %probability of influx
mySystem.params.M1influxRate = 1;   %number of cells per influx
mySystem.params.M1cellNo = 0;       %number of initial cells
mySystem.params.M1engagementDuration = 60; %number of steps the immune cell is engaged

%DEFINE MACROPHAGE TYPE 2 PARAMETERS
mySystem.params.M2kmax = 11;   %killing capacity 
mySystem.params.M2pkill = 0.0381; %probability of killing 
mySystem.params.M2pkill = 0.0381; %probability of mutated cell promotion
mySystem.params.M2pmig = 0.8001; %probability of migration
mySystem.params.M2pdeath = 0.0147; %probability of death
mySystem.params.M2rwalk = 0.8;      %random influence on movement
mySystem.params.M2speed = 40;       %speed of movement
mySystem.params.M2influxProb = 0; %probability of influx
mySystem.params.M2influxRate = 1;   %number of macrophages each influx
mySystem.params.M2cellNo = 0;       %number of initial cells
mySystem.params.M2TUadd = 0.2985; %Tumor promotion: added value to proliferation probability of tumor cells
mySystem.params.M2engagementDuration = 60; %number of steps immune cell is engaged
mySystem.params.M2emt = 0.75; %Chance of the M2 macrophage initating EMT in tumor cell

%DEFINE FIBROBLAST PARAMETERS
mySystem.params.Fpprol = 0; %probability of proliferation
mySystem.params.Fpmig = 0; %probability of migration
mySystem.params.Fpdeath = 0; %probability of death
mySystem.params.Fpmax = 4;  %initial proliferation capacity
mySystem.params.Frwalk = 0; %random influence on movement
mySystem.params.Fdiff = 0.75; %Chance of fibroblast turning into CAF when in range of mutated cell

% DEFINE CAF PARAMETERS
mySystem.params.CFpprol = (0.0838/10); %probability of proliferation
mySystem.params.CFpmig = 0.2; %probability of migration 
mySystem.params.CFpdeath = 0.0054; %Probability of death 
mySystem.params.CFpmax = 4; %Initial proliferation capacity 
mySystem.params.CFrwalk = 0.5; %random influence on movement 
mySystem.params.CFcellNo = 0; %Initial cell number
mySystem.params.CFemt = 0.75; %Chance of initiating EMT in tumor cells when CAF is close enough
mySystem.params.CFprom = 0.5*mySystem.params.TUpprol;  %factor with which tumorcell proliferation is multiplied by CAFs
mySystem.params.CFmmp = 0.9; %Chance of creating a path in stroma for mutated cells to move through (invade surroundings)
mySystem.params.CFmdiff = 0.5; %Probability of CAFs affecting macrophage differentiation

%DEFINE BASAL CELL PARAMETERS
mySystem.params.Cpprol =0.2; %proliferation probability 
mySystem.params.Cpmig = 0; %migration probability 
mySystem.params.Crwalk = 0; %random influence on movement
mySystem.params.Cpdeath = 0; %probability of death 
mySystem.params.Cpmax = 10; %initial proliferation capacity 

%DEFINE BASEMENT MEMBRANE
mySystem.params.Bpdeath = 0; %probability of death
mySystem.params.Bpprol = 0; %Proliferation probability 

%DEFINE DIMENSION PARAMETERS
%to create NxM grid
mySystem.grid.N=space(1);
mySystem.grid.M=space(2);

%DEFINE SOME CONSTANTS
cnst.nSteps = 800;                  %number of steps in the simulation
cnst.video = true;                 %Define whether or not you want to show the video 
cnst.newSystem = true;              %initialize new system or use previously defined
