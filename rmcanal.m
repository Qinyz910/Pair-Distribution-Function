
close all
clear

fileroot = '15k' ;
rmin = 0.5 ;
rmax = 2.0 ;
Deltar = 0.02 ;
Deltatheta = 1 ;

% Plot bond length distribution
filename = strcat(fileroot, '.rmcbonds') ;
data = dlmread(filename,'',0,5) ;
lengths = data ;
r = [rmin : Deltar : rmax] ;
y = hist(lengths,r) ;
figure(1)
plot(r,y,'ko','markerfacecolor','k')
xlabel('Distance (Ang)')
ylabel('Number')
goodplot

% Plot bond angle distribution
filename = strcat(fileroot, '.rmcangles') ;
data = dlmread(filename,'',0,7) ;
angles = data(:,3) ;
theta = [0 : Deltatheta : 180] ;
y = hist(angles,theta) ;
figure(2)
plot(theta,y,'ko','markerfacecolor','k')
xlabel('Angle (бу)')
ylabel('Number')
goodplot

print('-f1','bonds.eps','-depsc')
print('-f2','angles.eps','-depsc')

