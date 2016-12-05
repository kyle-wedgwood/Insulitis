% Script to plot some data

figure;
hold on;
set(gca,'Fontname','Times');
set(gca,'Fontsize',64);

j = 1;

mode = 'T cell';
cm   = [ 0       0.4470 0.7410 ;
         0.8500  0.3250 0.0980 ;
         0.9290  0.6940 0.1250 ];

if strcmp(mode,'beta cell');

  h = gobjects(1,2);
  
  while 1

    filename = input('Enter filename to load or done to finish - ','s');

    if strcmp(filename,'done');
      break;
    end

    load(sprintf('%s.mat',filename));

    T = linspace(0,size(betaCellMass,1),size(betaCellMass,1))';
    
    % Plot error bars
    sigma = std(betaCellMass(:,1:noTrials),0,2);
    X = [T;flipud(T)];
    Y = [betaCellMass(:,noTrials+1)+sigma;flipud(betaCellMass(:,noTrials+1)-sigma)];
    fill(X,Y,cm(mod(j-1,3)+1,:),'EdgeColor','none','FaceAlpha',0.3);
    h(j) = plot(T,betaCellMass(:,noTrials+1),'LineWidth',5,'Color',cm(mod(j-1,3)+1,:));
    j = j+1;

  end

  set(gca,'Xtick',0:2000:6000);
  set(gca,'ylim',[0,betaCellMass(1,noTrials+1)]);
  YTickLabels = 0:50:100;
  YTick = linspace(0,betaCellMass(1,noTrials+1),numel(YTickLabels));
  set(gca,'Ytick',YTick);
  set(gca,'Yticklabel',YTickLabels);
  
elseif strcmp(mode,'T cell')
  
  h = gobjects(1,2);
  
  while 1

    filename = input('Enter filename to load or done to finish - ','s');
    
    if strcmp(filename,'done');
      break;
    end

    load(sprintf('%s.mat',filename));
    
    betaCellMean = betaCellMass(:,noTrials+1);
    
    % Collect statistics
    sigmaT = zeros(size(betaCellMass,1),1);
    dataT  = zeros(size(betaCellMass,1),noTrials);
    sigmaB = zeros(size(betaCellMass,1),1);
    dataB  = zeros(size(betaCellMass,1),noTrials);
    
    for k = 1:noTrials
      [~,ind] = unique(betaCellMass(:,k));
      dataT(:,k) = interp1(betaCellMass(ind,k),noTInterior(ind,k),betaCellMean);
      dataB(:,k) = interp1(betaCellMass(ind,k),noBInterior(ind,k),betaCellMean);
    end
    for k = 1:size(betaCellMass,1)
      ind = ~isnan(dataT(k,:));
      sigmaT(k) = std(dataT(k,ind));
      sigmaB(k) = std(dataB(k,ind));
    end

    X = [betaCellMean;flipud(betaCellMean)];
    Y = [noTInterior(:,noTrials+1)-sigmaT;flipud(noTInterior(:,noTrials+1)+sigmaT)];
    fill(X,Y,cm(mod(j-1,3)+1,:),'EdgeColor','none','FaceAlpha',0.4);
    h(j,1) = plot(betaCellMass(:,noTrials+1),noTInterior(:,noTrials+1),'LineWidth',10,'Color',cm(mod(j-1,3)+1,:));
    j = j+1;
    
    Y = [noBInterior(:,noTrials+1)-sigmaB;flipud(noBInterior(:,noTrials+1)+sigmaB)];
    fill(X,Y,cm(mod(j-1,3)+1,:),'EdgeColor','none','FaceAlpha',0.4);
    j = 1;
    h(j,1) = plot(betaCellMass(:,noTrials+1),noTInterior(:,noTrials+1),'LineWidth',10,'Color',cm(mod(j-1,3)+1,:));
    j = j+1;
    h(j,2) = plot(betaCellMass(:,noTrials+1),noBInterior(:,noTrials+1),'LineWidth',10,'Color',cm(mod(j-1,3)+1,:));
    
  end
  
  set(gca,'xdir','reverse');
  set(gca,'xlim',[betaCellMass(end,noTrials+1),betaCellMass(1,noTrials+1)]);
  endPercentage = betaCellMass(end,noTrials+1)/betaCellMass(1,noTrials+1)*100;
  XTickLabels = 0:50:100;
  XTick = linspace(0,betaCellMass(1,noTrials+1),numel(XTickLabels));
  ylim = get(gca,'ylim');
  set(gca,'ylim',[0,ylim(2)]);
  npts = 4;
  YTickSpacing = floor(ylim(2)/npts);
  YTick = ylim(1):YTickSpacing:ylim(2);
  set(gca,'Xtick',XTick);
  set(gca,'Xticklabel',XTickLabels);
  set(gca,'YTick',YTick);
  
end