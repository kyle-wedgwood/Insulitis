% script to produce figures

close all;
filename = 'IBMTrialsCycleNoRepair';

load(sprintf('%s.mat',filename));

T = (0:parameters.nsteps/parameters.nPlot-1)'*100*parameters.dt;

figure(1);
[ax,h1,h2] = plotyy(T,noTInterior(:,end),T,betaCellMass(:,end));
hold on;
h3 = plot(T,noBInterior(:,end),'Parent',gca);

set(h1,'Linewidth',4);
set(h2,'color','k','Linewidth',4);
set(h3,'Linewidth',4);

xlim([0,400]);
set(gca,'XTick',0:100:400,'Fontname','Times-Roman','Fontsize',24);
set(ax(2),'Fontname','Times-Roman','Fontsize',24);

xlabel('Time (days)','Fontname','Times-Roman','Fontsize',32);
ylabel('Cell density','Fontname','Times-Roman','Fontsize',32);

saveFilename = sprintf('~/Dropbox/Insulitis/Figs/%sTimeCourse.eps',filename);
saveas(gcf,saveFilename,'epsc');

%%
figure(2);
h = plot(betaCellMass(:,end),noTInterior(:,end),betaCellMass(:,end),noBInterior(:,end));
set(gca,'XTick',0:100:300,'YTick',0:1:4,'Xdir','reverse','Fontname','Times-Roman','Fontsize',24);

for i = 1:length(h)
  set(h(i),'Linewidth',4);
end

xlabel('Beta cell density','Fontname','Times-Roman','Fontsize',32);
ylabel('Lymphocyte density','Fontname','Times-Roman','Fontsize',32);

saveFilename = sprintf('~/Dropbox/Insulitis/Figs/%sBeta.eps',filename);
saveas(gcf,saveFilename,'epsc');