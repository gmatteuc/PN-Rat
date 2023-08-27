function [hmu,hstd] = plot_shaded_mu_std(ax,xvals,mu,std,color,alpha)

hmu=plot(ax,xvals,mu,'-','Color',color,'Linewidth',3);
uE = mu+std;
lE = mu-std;
x = xvals;
yP = [lE,fliplr(uE)];
xP = [x,fliplr(x)];
idx2rm=isnan(yP);
xP(idx2rm)=[];
yP(idx2rm)=[];
hstd = patch(ax, xP,yP,1,'facecolor',color,...
    'edgecolor','none',...
    'facealpha',alpha);

end

