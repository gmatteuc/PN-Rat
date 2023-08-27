function [ handle] = plot_PGD_distribution( goopgdiff, target_group )

% [ handle] = plot_PGD_distribution( goopgdiff, target_group )
%
% plot PGD distribution
%-----------------------------------------------------------------------

binnum=20; %12

ff=figure;
set(ff,'Position',[10,10,1500,1000]);
[v1,c1]=hist(goopgdiff,binnum);
bar(c1,v1/sum(v1));
hold on
plot(c1,v1/sum(v1),'--b','LineWidth',3)
goopgdiff(isnan(goopgdiff))=0;
plot(median(goopgdiff).*ones(1,11),0:0.05:0.5,'--k','LineWidth',2.5)
ylimit=get(gca,'ylim');
xlimit=get(gca,'xlim');
te=text(-0.8*xlimit(2),0.40*ylimit(2),['median PGD = ',num2str(median(goopgdiff),'%.1f'),' Hz'],'FontSize',10);
xlim([-14,14])
ylim([0,0.45])
title(['pattern-grating difference distribution ',target_group])
set(gca, 'box','on');
axis square
hold off
fffname=['pgd distribution ',target_group];
set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf,fffname, 'jpg')


end

