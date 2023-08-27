function handle = plot_polar_mu_std(handle,angvec,muvec,stdvec,colortouse,alphatouse,linestyletouse,maxradiustouse)

% get inputs and prepare them for patch
thetavec=deg2rad([angvec,angvec(1)])';
rhovec=deg2rad([muvec,muvec(1)])';
rhovec_upper=rhovec+deg2rad([stdvec,stdvec(1)])';
rhovec_lower=rhovec-deg2rad([stdvec,stdvec(1)])';
tempang=[thetavec;thetavec];
temprad=[rhovec_upper;rhovec_lower];
[tempx,tempy] = pol2cart(tempang,temprad);
% draw patch and line
hold on;
patch(handle,'XData',tempx,'YData',tempy,...
    'FaceColor',colortouse,'FaceAlpha',alphatouse,...
    'EdgeColor',colortouse,'EdgeAlpha',0)
p1=polar(handle,thetavec,rhovec); %#ok<POLAR>
set(p1,'LineWidth',3); set(p1,'Color',colortouse);
set(p1,'LineStyle',linestyletouse);
% refine plot
axis equal
set(handle,'ycolor',[1,1,1])
set(handle,'xcolor',[1,1,1])
yticks(handle,'')
xticks(handle,'')
sizecircle=maxradiustouse; %max([ylimrange,xlimrange]);
pos = [-sizecircle/2 -sizecircle/2 sizecircle sizecircle];
rectangle(handle,'Position',pos,'Curvature',[1 1],'linewidth',2);
ylim(handle,[-sizecircle/2,sizecircle/2])
xlim(handle,[-sizecircle/2,sizecircle/2])
ylimused=get(gca,'ylim');
xlimused=get(gca,'xlim');
plot(handle,[0,0],ylimused,'--','linewidth',2,'color',[0.5,0.5,0.5])
plot(handle,xlimused,[0,0],'--','linewidth',2,'color',[0.5,0.5,0.5])

end

