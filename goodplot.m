function goodplot(fontsize,papersizex,papersizey, margin)
% function which produces a nice-looking plot
% and sets up the page for nice printing
if nargin == 0  % nargin is the number of incoming arguments
papersizex = 8;
papersizey = 6;
margin = 0.5;
fontsize = 18;
elseif nargin == 1
papersizex = 8;
papersizey = 6;
margin = 0.5;
elseif nargin == 2
papersizey = papersizex*0.75;
margin = 0.5;
elseif nargin == 3
margin = 0.5;
end
set(get(gca,'xlabel'),'FontName','Helvetica','FontSize', fontsize, 'FontWeight', 'Normal');
set(get(gca,'ylabel'),'FontName','Helvetica','FontSize', fontsize, 'FontWeight', 'Normal');
set(get(gca,'title'),'FontName','Helvetica','FontSize', fontsize, 'FontWeight', 'Normal');
set(gca,'TickLength',[0.02,0.0])
set(gca,'LineWidth',1);
set(gca,'FontSize',fontsize);
set(gca,'FontWeight','Normal');
set(gcf,'color','w');
set(gcf,'PaperPositionMode','Manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [papersizex papersizey]);
set(gcf,'PaperPosition',[margin margin papersizex-2*margin papersizey-2*margin]);
