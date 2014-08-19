function CombineTwoFig
close all;
fig1=open('LinearNor.fig');
fig2=open('SplineNor.fig');
ax1=get(fig1,'Children');
ax2=get(fig2,'Children');
for i=1:numel(ax2)
ax2Children=get(ax2(i),'Children');
copyobj(ax2Children,ax1(i));
end