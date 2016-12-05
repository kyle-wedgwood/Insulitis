% modelling forces acting on cells
function [F,G,H] = modelForcesAllCellTypesFaster(~,NT,NB,pos,betaCellPos,p,chemokine,cellType,cellRadii,theta,membrane,apoptotic,cellPars,interpFunc,activation,killingRate)

  debug = false;

  % unpack
  M = NT+NB;
  noBetaCells = size(betaCellPos,1)/2;
  F = zeros(2*M,1);
  iX = 1:2:2*M;
  iY = 2:2:2*M;
  x = pos(1:2:2*M);
  y = pos(2:2:2*M);
  betaX = betaCellPos(1:2:2*noBetaCells);
  betaY = betaCellPos(2:2:2*noBetaCells);

  isletRadius     = p(1);
  isletPosition   = p(2);

  sensingStrength = cellPars(:,1);
  viscosity       = cellPars(:,2);
  attractionStrength = cellPars(:,3);
  attractionRange    = cellPars(:,4);

  % decide which cells are internal and which are external
  r = sqrt((x-isletPosition).^2+y.^2);
  % index = [ r>isletRadius-cellRadii(1:NT+NB) , r<isletRadius+cellRadii(1:NT+NB) ]; % first column is exterior, second is interior

  % interaction with islet
  d = isletRadius-r;
  a = atan2(y,x-isletPosition);
  a(a<0) = 2*pi+a(a<0);

  membraneInterp = interpFunc(a(:,ones(size(theta)))-theta(ones(size(a)),:))*membrane';
  isletX = -100.0*heaviside(membraneInterp-0.1).*heaviside(1.0-abs(d)).*d.^5.*(x-isletPosition)./r;
  isletY = -100.0*heaviside(membraneInterp-0.1).*heaviside(1.0-abs(d)).*d.^5.*y./r;

  % local approximation to gradient of chemokine using centred differences
  dx = 0.01; % length of protrusions for gradient sensing
  dy = 0.01;
  R = cellRadii(NT+NB+1:end);
  chemoX = (chemokine(x+dx,y,R)-chemokine(x-dx,y,R))/(2*dx);
  chemoY = (chemokine(x,y+dy,R)-chemokine(x,y-dy,R))/(2*dy);

  %% interaction forces between cells
  %cellX = zeros(M,1);
  %cellY = zeros(M,1);
  %for i = 1:2
  %  cellCount = sum(index(:,i));
  %  localCellRadii = cellRadii(index(:,i));
  %  localCellType  = cellType(index(:,i));
  %  localX = x(index(:,i));
  %  localY = y(index(:,i));
  %  cellSumRadii = bsxfun(@plus,localCellRadii,localCellRadii');
  %  D = sqrt(bsxfun(@minus,localX,localX').^2+bsxfun(@minus,localY,localY').^2)-cellSumRadii;
  %  D(1:cellCount+1:cellCount*cellCount) = 1.0;
  %  attractionRange = attractionRange(:,ones(cellCount,1));
  %  cellX(index(:,i)) = cellX(index(:,i)) + sum(10.0*heaviside(-D).*D.^4.*bsxfun(@minus,xLocal,xLocal')./(D+cellSumRadii),2)...
  %      +2*attractionStrength.*sum(bsxfun(@ne,cellType,cellType')...
  %      .*bsxfun(@minus,xLocal,xLocal').*D.*exp(-D.^2./attractionRange.^2)./(attractionRange.^2*(D+cellSumRadii)),2);
  %  cellY(index(:,i)) = cellY(index(:,i)) + sum(10.0*heaviside(-D).*D.^4.*bsxfun(@minus,yLocal,yLocal')./(D+cellSumRadii),2)...
  %      +2*attractionStrength.*sum(bsxfun(@ne,cellType,cellType')...
  %      .*bsxfun(@minus,yLocal,yLocal').*D.*exp(-D.^2./attractionRange.^2)./(attractionRange.^2*(D+cellSumRadii)),2);
  %  H = sum(heaviside(0.5-D(1:NT,NT+1:NT+NB)),2)-0.5*activation;
  %end

  % interaction forces between cells
  cellSumRadii = bsxfun(@plus,cellRadii(1:NT+NB),cellRadii(1:NT+NB)');
  D = sqrt(bsxfun(@minus,x,x').^2+bsxfun(@minus,y,y').^2)-cellSumRadii;
  D(eye(M)==1) = 1.0;
  attractionRange = attractionRange(:,ones(M,1));
  cellX = sum(10.0*heaviside(-D).*D.^4.*bsxfun(@minus,x,x')./(D+cellSumRadii),2)...
      +2*attractionStrength.*sum(bsxfun(@ne,cellType,cellType').*bsxfun(@minus,x,x').*D.*exp(-D.^2./attractionRange.^2)./(attractionRange.^2*(D+cellSumRadii)),2);
  cellY = sum(10.0*heaviside(-D).*D.^4.*bsxfun(@minus,y,y')./(D+cellSumRadii),2)...
      +2*attractionStrength.*sum(bsxfun(@ne,cellType,cellType').*bsxfun(@minus,y,y').*D.*exp(-D.^2./attractionRange.^2)./(attractionRange.^2*(D+cellSumRadii)),2);
  if NT>0
    H = sum(heaviside(0.5-D(1:NT,NT+1:NT+NB)),2)-0.5*activation;
  else
    H = 0.0;
  end
  
  % interaction forces between immune cells and beta-cells
  betaCellX = zeros(M,1);
  betaCellY = zeros(M,1);
  interior  = r<isletRadius+cellRadii(1:NT+NB);
  interiorT = cellType(interior)==1;
  if any(interior)
    cellSumRadii = bsxfun(@plus,cellRadii(interior),cellRadii(NT+NB+1:end)');
    D = sqrt(bsxfun(@minus,x(interior),betaX').^2+bsxfun(@minus,y(interior),betaY').^2)-cellSumRadii;
    betaCellX(interior) = sum(10.0*heaviside(-D).*D.^4.*bsxfun(@minus,x(interior),betaX')./(D+cellSumRadii),2);
    betaCellY(interior) = sum(10.0*heaviside(-D).*D.^4.*bsxfun(@minus,y(interior),betaY')./(D+cellSumRadii),2);
    G = sum(bsxfun(@times,killingRate(interior(1:NT)),heaviside(0.5-D(interiorT,:))),1)'-0.5*apoptotic;
  else
    G = zeros(noBetaCells,1);
  end

  % resolve all forces
  F(iX) = (sensingStrength.*chemoX+isletX+cellX+betaCellX)./viscosity;
  F(iY) = (sensingStrength.*chemoY+isletY+cellY+betaCellY)./viscosity;

  if debug
    figure(2);
    ind = 1:size(pos,1)/2;
    plot(ind,sensingStrength*sqrt(chemoX.^2+chemoY.^2),ind,sqrt(isletX.^2+isletY.^2),ind,sqrt(cellX.^2+cellY.^2));
    drawnow;
  end

end
