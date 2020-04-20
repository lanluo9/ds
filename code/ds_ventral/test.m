t = 0 : .01 : 2*pi;
polarplot(t, sin(2*t).*cos(2*t),'Color',[0 1 0]);
hold on
polarplot(t, sin(2*t).*cos(2*t)+1,'Color',[.6 0 0]);
polarplot(t, sin(2*t).*cos(2*t)+2,'Color',[.6 0 0]);
% blue - red - orange
