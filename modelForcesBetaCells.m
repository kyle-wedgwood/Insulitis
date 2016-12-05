% modelling forces acting on beta-cells only
% V3 - now including beta cells as well
function F = modelForcesBetaCells(~,pos,p,cellRadii)

  debug = false;
  
  M = size(pos,1)/2;

  % unpack
  F = zeros(2*M,1);
  iX = 1:2:2*M;
  iY = 2:2:2*M;
  x = pos(1:2:2*M);
  y = pos(2:2:2*M);

  isletRadius     = p(1);
  isletPosition   = p(2);
  
  % interaction forces between cells
  cellSumRadii = bsxfun(@plus,cellRadii,cellRadii');
  D = sqrt(bsxfun(@minus,x,x').^2+bsxfun(@minus,y,y').^2)-cellSumRadii;
  D(eye(M)==1) = 1.0;
  cellX = sum(10.0*heaviside(-D).*D.^4.*bsxfun(@minus,x,x')./(D+cellSumRadii),2);
  cellY = sum(10.0*heaviside(-D).*D.^4.*bsxfun(@minus,y,y')./(D+cellSumRadii),2);

  % interaction with islet
  r = sqrt((x-isletPosition).^2+y.^2);
  d = r-isletRadius;
  isletX = -100.0*heaviside(1.2-abs(d)).*d.^4.*(x-isletPosition)./r;
  isletY = -100.0*heaviside(1.2-abs(d)).*d.^4.*y./r;

  % resolve all forces
  F(iX) = isletX+cellX;
  F(iY) = isletY+cellY;

end
