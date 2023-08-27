function [ handle1,handle2 ] = plot_pattern_component( Zc, Zp, PI, Zcgen, Zpgen, PIgen, goodprds,selectedsi, target_group )

% [ handle1,handle2 ] = plot_pattern_component( Zc, Zp, PI,goodprds,selectedsi, target_group )
%
% plot pattern and component analysis
%-----------------------------------------------------------------------

nnV=cell2mat(selectedsi);
nnV=nnV(:,1);

% ------- Pattern vs. Component plot -------

handle1=figure;
set(handle1,'Position',[10,10,1500,1000]);
par=winter;
% plot every neuron
plot(Zc,Zp,'.','MarkerSize',45,'Color',par(1,:));
hold on
% plot only DS and plaid responsive ones
for nnn=1:length(goodprds)
    if goodprds(nnn)
        plot(Zc(nnn),Zp(nnn),'.','MarkerSize',25,'Color',par(end,:));
        ylimit=get(gca,'ylim');
        te=text(Zc(nnn),Zp(nnn)+0.07*diff(ylimit),['n = ',num2str(nnV(nnn))],'FontSize',15);
    else
    end
end
line([0 5], [1.28 6.28],'LineWidth',1.5,'Color','k');
line([1.28 6.28], [0 5],'LineWidth',1.5,'Color','k');
line([1.28 1.28], [-4 0],'LineWidth',1.5,'Color','k');
line([-4 0], [1.28 1.28],'LineWidth',1.5,'Color','k');
xlabel('Zc'); ylabel('Zp'); title('Pattern vs. Component');
legend([target_group,' neurons'],'selected neurons','Location','NorthEast')
% save
saveas(handle1,[target_group,'_pattern_component_Z.eps'], 'eps');
set(handle1, 'PaperPositionMode', 'auto'); saveas(handle1,[target_group,'_pattern_component_Z.jpg'], 'jpg');

% ------- Pattern vs. Component plot -------

handle1=figure;
set(handle1,'Position',[10,10,1500,1000]);
par=winter;
% plot every neuron
plot(Zcgen,Zpgen,'.','MarkerSize',45,'Color',par(1,:));
hold on
% plot only DS and plaid responsive ones
for nnn=1:length(goodprds)
    if goodprds(nnn)
        plot(Zcgen(nnn),Zpgen(nnn),'.','MarkerSize',25,'Color',par(end,:));
        ylimit=get(gca,'ylim');
        te=text(Zcgen(nnn),Zpgen(nnn)+0.07*diff(ylimit),['n = ',num2str(nnV(nnn))],'FontSize',15);
    else
    end
end
line([0 5], [1.28 6.28],'LineWidth',1.5,'Color','k');
line([1.28 6.28], [0 5],'LineWidth',1.5,'Color','k');
line([1.28 1.28], [-4 0],'LineWidth',1.5,'Color','k');
line([-4 0], [1.28 1.28],'LineWidth',1.5,'Color','k');
xlabel('Zcgen'); ylabel('Zpgen'); title('Pattern vs. Component (generalized)');
legend([target_group,' neurons'],'selected neurons','Location','NorthEast')
% save
saveas(handle1,[target_group,'_generalized_pattern_component_Z.eps'], 'eps');
set(handle1, 'PaperPositionMode', 'auto'); saveas(handle1,[target_group,'_generalized_pattern_component_Z.jpg'], 'jpg');

% ------------ PI distribution ------------

% plot for every neuron
handle2=figure; [n1,x1] = hist(PI);
p1 = bar(x1,n1/sum(n1));
xlim([-8,8]); ylim([0,0.3]);
xlabel('PI');
ylabel('fraction of neurons'); title([target_group,' PI distribution - n=',num2str(numel(nnV))]);
% save
set(gcf, 'PaperPositionMode', 'auto'); saveas(handle2, [target_group,' pattern_index.jpg'], 'jpg');

% plot for every neuron
handle2=figure; [n1,x1] = hist(PIgen);
p1 = bar(x1,n1/sum(n1));
xlim([-8,8]); ylim([0,0.3]);
xlabel('PI generalized');
ylabel('fraction of neurons'); title([target_group,' PI generalized distribution - n=',num2str(numel(nnV))]);
% save
set(gcf, 'PaperPositionMode', 'auto'); saveas(handle2, [target_group,' generalized_pattern_index.jpg'], 'jpg');


end

