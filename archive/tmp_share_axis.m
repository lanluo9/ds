% ax = tight_subplot(3,1,0,0.1,.1);hold on
% set(ax,'color','none')
% % set(ax(3:4),'xcolor','none')
% % set(ax([1 3]),'ycolor','none')
% % set([ax.XRuler],'tickdirection','both')
% % set([ax.YRuler],'tickdirection','both')
% for i=1:3
%     axes(ax(i));hold on
%     scatter(rand(10,1),rand(10,1))
% end
% linkaxes(ax,'xy')

ax = tight_subplot(11,1,0,0.1,0.1); hold on
set(ax,'color','none')
% set(ax(3:4),'xcolor','none')
% set(ax([1 3]),'ycolor','none')
% set([ax.XRuler],'tickdirection','both')
% set([ax.YRuler],'tickdirection','both')
for i=1:2
    axes(ax(i)); hold on
    axes('Position', [0.1 + 0.4*i, 0.1, 0.4, 0.8]);
%     axes('Position', [0.5, 0.1, 0.4, 0.8]); 
    scatter(rand(10,1),rand(10,1))
end
linkaxes(ax,'y')