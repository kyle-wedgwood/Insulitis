% wrapper to collect statistics for the IBM model

% T-cells
parameters.T.id = 1;
parameters.T.N  = 30;
parameters.T.cellRadius      = 1.2;
parameters.T.sensingStrength = 600.0;
parameters.T.viscosity       = 1.0;
parameters.T.attractionStrength = 0.0;
parameters.T.attractionRange    = 1.0;
parameters.T.activationRate     = 0.2;
parameters.T.activatedSensitivityBoost = 400.0;
parameters.T.killingRate = 0.2; 
parameters.T.activatedKillingRateBoost = 0.9;
parameters.T.cellCycleLength   = 56000;
parameters.T.activatedCellCycleDecrease = 1;

% B-cells
parameters.B.id = 2;
parameters.B.N  = 30;
parameters.B.cellRadius      = 1.2;
parameters.B.sensingStrength = 700;
parameters.B.viscosity       = 1.0;
parameters.B.attractionStrength = 1.0;
parameters.B.attractionRange    = 1.0;
parameters.B.cellCycleLength = 56000;

% beta-cells
parameters.beta.id = 3;
parameters.beta.N  = 30;
parameters.beta.cellRadius    = 1.6;
parameters.beta.apoptosisRate = 0.1;

% islet
parameters.islet.Radius   = 30.0;
parameters.islet.Position = 0.0;

% domain
parameters.domain.Height = 80.0;
parameters.domain.Width  = 250.0;
parameters.domain.Radius = 100.0;

% noise
parameters.sigma = 1.0;
parameters.sigmaSensing = 0.5;

% membrane
parameters.repairRate       = 0.01;
parameters.degradationRate  = 0.1;
parameters.degradationRange = 1;

% solver options
parameters.tstart = 0.0;
parameters.tfinal = 600.0;
parameters.dt     = 0.001;
parameters.nsteps = ceil((parameters.tfinal-parameters.tstart)/parameters.dt);

% chemokine fun
parameters.a1 = 50;
parameters.a2 = 2;
parameters.chemoStrength1 = 0.5;
parameters.chemoStrength2 = 10.0;

% output
parameters.nPlot = 100;
parameters.folderName = 'Test';

noTrials = 1;
nSave = parameters.nsteps/parameters.nPlot;
betaCellMass = zeros(nSave,noTrials);
noTInterior = zeros(nSave,noTrials);
noBInterior = zeros(nSave,noTrials);
noMembraneHoles = zeros(nSave,noTrials+1);
noActivatedTCells = zeros(nSave,noTrials+1);

for i = 1:noTrials
  tic;
  data = IBM(parameters);
  toc;
  betaCellMass(:,i) = data(:,1);
  noTInterior(:,i)  = data(:,2);
  noBInterior(:,i)  = data(:,3);
  noMembraneHoles(:,i) = data(:,4);
  noActivatedTCells(:,i) = data(:,5);

  fprintf('\nSimulation %d of %d complete\n',i,noTrials);
end

clear nSave data

% do averaging
betaCellMass(:,noTrials+1) = mean(betaCellMass(:,1:noTrials),2);
noTInterior(:,noTrials+1)  = mean(noTInterior(:,1:noTrials),2);
noBInterior(:,noTrials+1)  = mean(noBInterior(:,1:noTrials),2);
noMembraneHoles(:,noTrials+1) = mean(noMembraneHoles(:,1:noTrials),2);
noActivatedTCells(:,noTrials+1) = mean(noActivatedTCells(:,1:noTrials),2);
filename = 'IBMTrials.mat';
if exist(filename,'file')~=2
  save(filename);
else
  i = 1;
  while 1
    altFilename = sprintf('IBMTrials%d.mat',i);
    if exist(altFilename,'file')~=2
      save(altFilename);
      break;
    end
    i = i+1;
  end
end