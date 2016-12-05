% plot membrane with varying color
function plotMembrane(theta,isletPosition,isletRadius,membrane)

  % hack to get BW colors to work
  theta = [theta,0];
  
  x = isletPosition+isletRadius*cos(theta);
  y = isletRadius*sin(theta);
  
  membrane = [heaviside(membrane-0.1),0.0];
    
  surface([x;x],[y;y],zeros(2,length(x)),zeros(2,length(x)),...
        'facecol','no',...
        'edgecol','interp',...
        'AlphaData',[membrane;membrane],...
        'edgeAlpha','interp',...
        'linew',2);
    
end