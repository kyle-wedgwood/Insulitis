% Individual based model for insulitis first example
% V13 - adding in cell turnover
%     - changing domain shape
%     - making cells appear anywhere in the extra islet space

function data = IBMV13(parameters)

  T = parameters.T;
  B = parameters.B;
  beta = parameters.beta;

  % islet
  isletRadius   = parameters.islet.Radius;
  isletPosition = parameters.islet.Position;

  % domain
  domainRadius = parameters.domain.Radius;

  % noise
  sigma = parameters.sigma;

  % initialise membrane
  theta = 0.0:0.01:2*pi;
  membrane = ones(size(theta));
  repairRate       = parameters.repairRate;
  degradationRate  = parameters.degradationRate;
  degradationRange = parameters.degradationRange;

  % activation of T-cells
  activationRate = parameters.T.activationRate;

  % killing rate
  killingRate = parameters.T.killingRate*ones(T.N,1);

  % solver options
  tstart = parameters.tstart;
  tfinal = parameters.tfinal;
  dt     = parameters.dt;
  nsteps = parameters.nsteps;

  % pack parameters
  p = [ isletRadius ;
      isletPosition ;
      domainRadius ];

  % store output
  trajectories = zeros(2*T.N+2*B.N,nsteps);
  nPlot = parameters.nPlot;
  nSave = nsteps/nPlot;
  data  = zeros(nSave,5);

  % auxiliary variables
  iX = 1:2:2*(T.N+B.N);
  iY = 2:2:2*(T.N+B.N);

  % radii sum
  distfun = @(x,y) x+y;

  % initialise system
  while 1
    x0 = rand(2*(T.N+B.N),1);
    cellRadii = [ T.cellRadius*ones(T.N,1); B.cellRadius*ones(B.N,1) ];
    cellType  = [ T.id*ones(T.N,1) ; B.id*ones(B.N,1) ];
    rho = (domainRadius-isletRadius-cellRadii(1:T.N+B.N)).*x0(1:2:end)+isletRadius+cellRadii(1:T.N+B.N);
    th  = 2*pi*x0(2:2:end);
    x0(1:2:end) = rho.*cos(th);
    x0(2:2:end) = rho.*sin(th);
    D = pdist([x0(1:2:end),x0(2:2:end)],'euclidean');
    if all(D>pdist(cellRadii,distfun))
      break;
    end
  end
%   cellCycleLength = [ones(T.N,1)*T.cellCycleLength;ones(B.N,1)*B.cellCycleLength]; % lifespan of cells
  cellCycleLength = [ T.cellCycleLength ; B.cellCycleLength ];
  counter = [randi(T.cellCycleLength,T.N,1);randi(B.cellCycleLength,B.N,1)]; % counter telling you how far along life cycle you are
  disp('Initial cell position configured.');

  % initialise islet
  noLayers = (isletRadius-beta.cellRadius)/(2*beta.cellRadius);
  dr = (isletRadius-beta.cellRadius)/noLayers;
  circum = 2*pi*(0:noLayers-1)'*dr;
  noPerLayer = floor((circum-2*beta.cellRadius)/(2*beta.cellRadius));
  dtheta = (circum-2*beta.cellRadius)./noPerLayer;
  beta.N = sum(noPerLayer(2:end))+1;
  y0 = zeros(2*beta.N,1);
  y0(1) = isletPosition;
  y0(2) = eps;
  j  = 1;
  for i = 2:noLayers
    xPos = repmat((i-1)*dr,[noPerLayer(i),1]).*cos(2*pi*(0:noPerLayer(i)-1)'/noPerLayer(i));
    yPos = repmat((i-1)*dr,[noPerLayer(i),1]).*sin(2*pi*(0:noPerLayer(i)-1)'/noPerLayer(i));
    y0(2*j+1:2:2*(j+noPerLayer(i))) = xPos+isletPosition;
    y0(2*j+2:2:2*(j+noPerLayer(i))) = yPos;
    j = j+noPerLayer(i);
  end
  cellRadii = [ cellRadii ; beta.cellRadius*ones(beta.N,1) ];
  shakeSteps = 5000;
  for i = 1:shakeSteps
    y0 = y0 + dt*modelForcesBetaCells(0,y0,p,cellRadii(T.N+B.N+1:end)) + sigma*sqrt(dt)*randn(2*beta.N,1);
  end
  disp('Islet configured.');

  % pack auxiliary parameters
  sensingStrength = [ T.sensingStrength ; B.sensingStrength ];
  viscosity       = [ T.viscosity ; B.viscosity ];
  attractionStrength = [ T.attractionStrength ; B.attractionStrength ];
  attractionRange = [ T.attractionRange ; B.attractionRange ];
  cellPars = [ sensingStrength(cellType) ,...
               viscosity(cellType) ,...
               attractionStrength(cellType) ,...
               attractionRange(cellType) ];

  % put randomness in sensing strength
  sigmaSensing = parameters.sigmaSensing;
  cellPars(:,1) = (1+sigmaSensing^2*randn(size(cellPars,1),1)).*cellPars(:,1);

  x0(iX) = x0(iX);
  trajectories(:,1) = x0;

  % distribute external beta cells
  noExternalBetaCells = 6;
  z0 = zeros(noExternalBetaCells,2);
  while 1
    rhoExternal   = 6*isletRadius + 4*isletRadius*rand(noExternalBetaCells,1);
    thetaExternal = 2*pi*rand(noExternalBetaCells,1);
    z0(:,1) = rhoExternal.*cos(thetaExternal);
    z0(:,2) = rhoExternal.*sin(thetaExternal);
    D = pdist(z0,'euclidean');
    if all(D>2*isletRadius);
      break;
    end
  end
  disp('External islets positioned.');

  % create chemokine function
  a1 = parameters.a1;
  a2 = parameters.a2;
  chemoStrength1 = parameters.chemoStrength1;
  chemoStrength2 = parameters.chemoStrength2;
  chemokine = @(x,y,R) chemoStrength1/beta.N*sum(bsxfun(@times,exp(-(bsxfun(@minus,x,y0(1:2:end)').^2+bsxfun(@minus,y,y0(2:2:end)').^2)/(2*a1^2)),R'),2)...;
                   + chemoStrength2/beta.N*sum(bsxfun(@times,exp(-(bsxfun(@minus,x,y0(1:2:end)').^2+bsxfun(@minus,y,y0(2:2:end)').^2)/(2*a2^2)),R'),2)...
                   + 0.5*chemoStrength1/(noExternalBetaCells*beta.N)*sum(exp(-(bsxfun(@minus,x,z0(:,1)').^2+bsxfun(@minus,y,z0(:,2)').^2)/(2*a1^2)),2); % external islets

  % initialise plot with circular domain
  npts = 200;
  xChemo = linspace(-domainRadius,domainRadius,npts);
  yChemo = linspace(-domainRadius,domainRadius,npts);
  [xC,yC] = meshgrid(xChemo,yChemo);

  % function for interpolation
  h = 2*pi/length(theta);
  interpFunc = @(x) sin(pi*x/h) ./( tan(x/2) .* (2*pi/h) );

  figure(1);
  hold on;
  maxChemoColor = max(chemokine(xC(:),yC(:),cellRadii(T.N+B.N+1:end))); % for caxis
  imagesc(xChemo,yChemo,1-reshape(chemokine(xC(:),yC(:),cellRadii(T.N+B.N+1:end)),[npts,npts])/maxChemoColor);colormap('gray');
  caxis([0,1]);
  plotMembrane(theta,isletPosition,isletRadius,membrane);
  plotCellsMultipleCellTypes([x0;y0],[cellType;beta.id*ones(beta.N,1)],cellRadii,1);
  axis equal;
  axis off;

  % apoptosis
  apoptotic             = zeros(beta.N,1);
  apoptoticCells        = zeros(beta.N,1);
  nonApoptoticCellIndex = (1:beta.N)';

  % T-cell activation
  activation     = zeros(T.N,1);
  activatedCells = zeros(T.N,1);

  % step & plot counter
  k = 0;
  imageBlockSize = 8000;
  disp('Simulation starting.');

  % for displaying purposes
  reverseStr = '';

  for i = 2:nsteps

    eta = randn(2*(T.N+B.N),1);
    [F,G,H] = modelForcesAllCellTypesFaster(0,T.N,B.N,x0,y0,p,chemokine,cellType,cellRadii,theta,membrane,apoptotic,cellPars,interpFunc,activation,killingRate);
    x1  = x0 + F*dt + sqrt(dt)*sigma*eta;
    if any(isnan(x1));
      pause;
    end
    xPos = x1(iX);
    yPos = x1(iY);

    % cells close to membrane
    D   = bsxfun(@minus,sqrt(bsxfun(@minus,xPos(cellType==1),isletPosition+isletRadius*cos(theta)).^2+...
          bsxfun(@minus,yPos(cellType==1),isletRadius*sin(theta)).^2),cellRadii(cellType==1));
    trajectories(:,i) = x1;
    membrane = (membrane+dt*repairRate)./(1+dt*repairRate+dt*degradationRate.*sum(exp(-D/degradationRange),1));
  %   membrane = membrane + dt*(repairRate*(1-membrane)-membrane*degradationRate.*sum(exp(-D/degradationRange),1));

    % apoptosis of beta cells
    apoptotic = apoptotic + dt*G;
    apoptoticCells = (apoptoticCells==1)|(apoptotic>1);
    cellRadii(find(apoptoticCells)+T.N+B.N) = 0.0;

    % activation of T-cells
    activation = activation + dt*activationRate*H;
    ind = (activation>1)&~activatedCells;
    cellPars(ind,1) = cellPars(ind,1) + T.activatedSensitivityBoost;
    killingRate(ind) = killingRate(ind) + T.activatedKillingRateBoost;
    activatedCells   = (activatedCells==1)|(activation>1);

    % immune cell death
    indCellDeath = find(counter>cellCycleLength(cellType));
    for j = 1:numel(indCellDeath)
      currentCellInd = indCellDeath(j);
      counter(currentCellInd) = 0;
      cellPars(currentCellInd,1) = (1+sigmaSensing^2*randn)*sensingStrength(cellType(currentCellInd));
      if currentCellInd<=T.N
        activation(currentCellInd) = 0.0;
        activatedCells(currentCellInd) = 0;
        killingRate(currentCellInd) = parameters.T.killingRate;
      end
      while 1
        th  = 2*pi*rand;
        rho = (domainRadius-isletRadius-cellRadii(currentCellInd))*rand+isletRadius+cellRadii(currentCellInd);
        x1(2*currentCellInd-1) = rho*cos(th);
        x1(2*currentCellInd) = rho*sin(th);
        if sum(sqrt((x1(2*currentCellInd-1)-x1(1:2:end)).^2 + (x1(2*currentCellInd)-x1(2:2:end)).^2) < cellRadii(currentCellInd) + cellRadii(1:T.N+B.N)) == 1
          break;
        end
      end
    end
    counter(1:T.N) = counter(1:T.N) + 1 + activatedCells*T.activatedCellCycleDecrease;
    counter(T.N+1:T.N+B.N) = counter(T.N+1:T.N+B.N)+1;

    % plot stuff
    if mod(i,nPlot)==0

      interior = sqrt((x1(1:2:end)-isletPosition).^2+x1(2:2:end).^2)<isletRadius;
      betaCellMass = sum(cellRadii(T.N+B.N+1:T.N+B.N+beta.N));
      noTInterior  = sum(interior(1:T.N));
      noBInterior  = sum(interior(T.N+1:T.N+B.N));
      noMembraneHoles = sum(membrane<0.1);
      data(k+1,:) = [betaCellMass,noTInterior,noBInterior,noMembraneHoles,sum(activatedCells)];

      if mod(k,imageBlockSize)==0
        if floor(k/imageBlockSize)==999
          break;
        end
        folderName = sprintf('Anim/%s/Block%03d',parameters.folderName,floor(k/imageBlockSize));
        if exist(folderName,'dir')~=7
          mkdir(folderName);
        end
      end
    
      figure(1);
      cla;
      imagesc(xChemo,yChemo,1-reshape(chemokine(xC(:),yC(:),cellRadii(T.N+B.N+1:end)),[npts,npts])/maxChemoColor);colormap('gray');
      plotMembrane(theta,isletPosition,isletRadius,membrane);
      plotCellsMultipleCellTypesActivated([x1;y0],[cellType;beta.id*ones(beta.N,1)],cellRadii,activatedCells,1);
      caxis([0,1]);
      axis equal;
      axis off;
      drawnow;

%     figure(2);
%     plot(theta,membrane);
% 
%     figure(3);
%     plot((1:beta.N),apoptotic,(1:beta.N),apoptoticCells);

     filename = sprintf('%s/image%07d.png',folderName,k);
     saveas(1,filename,'png');

     k = k+1;

     % display
     msg = sprintf('Done %d of %d', k, nSave);
     fprintf([reverseStr,msg]);
     reverseStr = repmat(sprintf('\b'), 1, length(msg));

    end

    % prepare for next step
    x0 = x1;

  end

end
