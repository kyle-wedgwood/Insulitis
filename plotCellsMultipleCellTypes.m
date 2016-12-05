% plot cells on figure handle
function plotCellsMultipleCellTypes(x,cellType,cellRadii,parentHandle)

  M  = size(x,1);
  cm = {'r','g','c'};
  theta = (0:0.01:2*pi)';
  xPos = x(1:2:M)';
  yPos = x(2:2:M)';

  for i = 1:3
    xdata = bsxfun(@plus,bsxfun(@times,cos(theta),cellRadii(cellType==i)'),xPos(cellType==i));
    ydata = bsxfun(@plus,bsxfun(@times,sin(theta),cellRadii(cellType==i)'),yPos(cellType==i));
    figure(parentHandle);
    plot(xdata,ydata,'Color',cm{i});
  end

end
