function [mySystem, cnst]=getSystemParams_variable(space,timestep)
%Sets the parameters needed for modeling the ABM (using a variable timestep, defined by the user)

mySystem.params.initialSeed = 1;   % initial random seed, default 1
mySystem.params.stromaPerm = 0.001; % stroma permeability in general (chance that an agent can take up a space that was stroma initially)
mySystem.params.stromaImPerm = 0.8; %stroma permeability for immune cells 
mySystem.params.Shape = 1; %Setting the shape of the acinus: 2 indicating ellips, 1 indicating circle 
mySystem.params.Number = 1; %Indicating the number of acini that you want 
mySystem.params.Fradius = 0; %Radius of the acinus
mySystem.params.Fellips= 0;  %Second 'radius' of the acinus, if you are plotting an ellips
mySystem.params.timestep = timestep; %The defined timestep by the user

%DEFINE TUMOR CELL PARAMETERS
mySystem.params.TUpprol = 0.02545*timestep; %probability of proliferation
mySystem.params.TUpmig = 0.029175*timestep; %probability of migration
mySystem.params.TUpdeath = 0.00071*timestep; %probability of death
mySystem.params.TUrwalk = 0.5;    %random influence on movement
mySystem.params.TUpmax = 4;   %initial proliferation capacity
mySystem.params.TUps = 0.024725*timestep; %probability of symmetric division
mySystem.params.TUcellNo = 1500;      %Number of initial luminal cells, only used when randomly seeding luminal cells (otherwise it plots the shape of the acinus)
mySystem.params.TUpmut = 0.000125*timestep; %probability of acquiring a mutation  
mySystem.params.TUantig = 0.005; %Added antigenicity per mutation 
mySystem.params.TUmutmax = 1; %maximum antigenicity a cell can have (i.e. TUmutmax/Tuantig = amount of mutations one cell can acquire) 
mySystem.params.TUwidth = 2; %Width of the layer of initial healthy luminal cells
mySystem.params.TUpkill = 0.025*timestep; %Chance of breaking down Basement membrane upon reaching the amount of required mutations
mySystem.params.TUthrshBM = (mySystem.params.TUantig*3); %Threshold of accumulated mutations value for breaking down BM
mySystem.params.TUpMdiff = 0.0667*timestep; %Chance that a tumor cell affects macrophage differentiation (via cytokines)
mySystem.params.TUthrshM = (mySystem.params.TUantig*4); %Needed threshold of accumulated mutations for differentiation effect
mySystem.params.TUpMMP = 0*timestep; %Chance that tumor cells create a path in the stroma to migrate to/through
mySystem.params.TUthrshMMP = (mySystem.params.TUantig*4); %Required amount of mutations before being able to break down stroma

%DEFINE MACROPHAGE TYPE 1 PARAMETERS
mySystem.params.M1kmax = 11;        %killing capacity 
mySystem.params.M1pkill = 0.00765*timestep; %probability of killing 
mySystem.params.M1pmig = 0.066675*timestep; %probability of migration
mySystem.params.M1pdeath = 0.001225*timestep; %probability of death
mySystem.params.M1rwalk = 0.8;      %random influence on movement
mySystem.params.M1speed = round(40/12)*timestep;       %speed of movement (amount of times a cell can move per timestep)
mySystem.params.M1influxProb = 0.025*timestep; %probability of influx
mySystem.params.M1influxRate = 1;   %number of cells per influx
mySystem.params.M1cellNo = 0;       %number of initial cells
mySystem.params.M1engagementDuration = 1.5*mySystem.params.M1speed; %number of steps immune cell is engaged

%DEFINE MACROPHAGE TYPE 2 PARAMETERS
mySystem.params.M2kmax = 11;   %killing capacity 
mySystem.params.M2pkill = 0.003175*timestep; %probability of killing 
mySystem.params.M2pkill = 0.003175*timestep; %probability of cell mutated promotion
mySystem.params.M2pmig = 0.066675*timestep;% %probability of migration
mySystem.params.M2pdeath = 0.001225*timestep; %probability of death
mySystem.params.M2rwalk = 0.8;      %random influence on movement
mySystem.params.M2speed = round(40/12)*timestep;       %speed of movement
mySystem.params.M2influxProb = 0; %probability of influx
mySystem.params.M2influxRate = 1;   %number of macrophages each influx
mySystem.params.M2cellNo = 0;       %number of initial cells
mySystem.params.M2TUadd = 0.024875*timestep; %Tumor promotion: added value to proliferation probability of tumor cells
mySystem.params.M2engagementDuration = 1.5*mySystem.params.M2speed; %number of steps immune cell is engaged
mySystem.params.M2emt = 0.0625*timestep; %Chance of the M2 macrophage initating EMT in tumor cell

%DEFINE FIBROBLAST PARAMETERS
mySystem.params.Fpprol = 0; %probability of proliferation
mySystem.params.Fpmig = 0; %probability of migration
mySystem.params.Fpdeath = 0; %probability of death
mySystem.params.Fpmax = 4; %initial proliferation capacity
mySystem.params.Frwalk = 0; %random influence on movement
mySystem.params.Fdiff = 0.0625*timestep; %Chance of fibroblast turning into CAF when in range of mutated cell

% DEFINE CAF PARAMETERS
mySystem.params.CFpprol = 0.00069833*timestep; %probability of proliferation
mySystem.params.CFpmig = 0.0166*timestep; %probability of migration
mySystem.params.CFpdeath = 0.00045*timestep; %Probability of death 
mySystem.params.CFpmax = 4; %Initial proliferation capacity 
mySystem.params.CFrwalk = 0.5; %random influence on movement 
mySystem.params.CFcellNo = 0; %Initial cell number
mySystem.params.CFemt = 0.0625*timestep; %Chance of initiating EMT in tumor cells when CAF is close enough
mySystem.params.CFprom = 0.5*mySystem.params.TUpprol;  %factor with which tumorcell proliferation is multiplied by CAFs
mySystem.params.CFmmp = 0.0625*timestep; %Chance of creating a path in stroma for mutated cells to move through (invade surroundings)
mySystem.params.CFmdiff = 0.0417*timestep; %Probability of CAFs affecting macrophage differentiation

%DEFINE BASAL CELL PARAMETERS
mySystem.params.Cpprol =0.05*timestep; %proliferation probability 
mySystem.params.Cpmig = 0; %migration probability 
mySystem.params.Crwalk = 0; %random influence on movement
mySystem.params.Cpdeath = 0.0025*timestep; %probability of death 
mySystem.params.Cpmax = 10; %initial proliferation capacity 

%DEFINE BASEMENT MEMBRANE
mySystem.params.Bpdeath = 0; %probability of death
mySystem.params.Bpprol = 0; %Proliferation probability 

%DEFINE DIMENSION PARAMETERS
%to create NxM grid
mySystem.grid.N=space(1);
mySystem.grid.M=space(2);

%DEFINE SOME CONSTANTS
cnst.nSteps = 100;                    %number of steps in the simulation
cnst.video = false;                   %Define whether or not you want to show the video 
cnst.newSystem = true;                %initialize new system or use previously defined
