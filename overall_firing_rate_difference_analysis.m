clear
close all
clc

% Still 2019 version... adapt to 2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set pars

pars = set_pars_PN();

listSessions = pars.listSessions;

code_folder=pars.code_folder;
addpath(code_folder);
data_folder=pars.processed_data_folder;
addpath(data_folder);

target_areas={'V1','LM','RL'};

outdirname=[data_folder,filesep,'2019_results',filesep,'friring_rate_difference_analysis_results'];
mkdir(outdirname)
oldfold=cd(outdirname);

%% load neuron selection results

% load indexing
load('Indexing.mat')
load('Tuning.mat')

% load each area's data
for area_idx=1:length(target_areas)
    
    % set area
    are=target_areas{area_idx};
    % load selection results
    sel_res_file=[are,'_selected_neurons',filesep,'selected_population_',are,'.mat'];
    results.(are)=load(sel_res_file);
    
end

for area_idx=1:length(target_areas)
    
    % set area
    are=target_areas{area_idx};
    
    % get current number of neurons
    current_n_neu=length(results.(are).selectedsi);
    
    % initialize max fr storage vector
    max_best_plaid_tc=zeros(current_n_neu,1);
    max_best_grating_tc=zeros(current_n_neu,1);
    max_best_tc=zeros(current_n_neu,1);
    
    % get max fr for gratings and plaids
    for neu_idx=1:current_n_neu
        k=results.(are).selectedsi{neu_idx}(1);
        i=results.(are).selectedsi{neu_idx}(2);
        j=results.(are).selectedsi{neu_idx}(3);
%         i=2; %SF 0.04
%         j=1; %TF 2
        best_plaid_tc=tuning_curve(:,k,i,j,2);
        best_grating_tc=tuning_curve(:,k,i,j,1);
        max_best_plaid_tc(neu_idx)=max(best_plaid_tc);
        max_best_grating_tc(neu_idx)=max(best_grating_tc);
        max_best_tc(neu_idx)=max(max_best_plaid_tc(neu_idx),max_best_grating_tc(neu_idx));
    end
    
    % identify cell types
    pattern_ds=((results.(are).goomlabel.*results.(are).goodprds)==2);
    component_ds=((results.(are).goomlabel.*results.(are).goodprds)==1);
    pattern_db=((results.(are).goomlabel.*results.(are).goodprdb)==2);
    component_db=((results.(are).goomlabel.*results.(are).goodprdb)==1);
    pattern_nds=((results.(are).goomlabel.*results.(are).goodpr)==2);
    component_nds=((results.(are).goomlabel.*results.(are).goodpr)==1);
    pattern_os=((results.(are).goomlabel.*results.(are).goodpros)==2);
    component_os=((results.(are).goomlabel.*results.(are).goodpros)==1);
    pattern_ob=((results.(are).goomlabel.*results.(are).goodprob)==2);
    component_ob=((results.(are).goomlabel.*results.(are).goodprob)==1);
    % get pattern index
    pi_ds=results.(are).goopi(find(results.(are).goodprds));
    pi_db=results.(are).goopi(find(results.(are).goodprdb));
    pi_nds=results.(are).goopi(find(results.(are).goodpr));
    pi_os=results.(are).goopi(find(results.(are).goodpros));
    pi_ob=results.(are).goopi(find(results.(are).goodprob));
    % get precomputed firing rate difference
    frdiff_ds=results.(are).goopgdiff(find(results.(are).goodprds));
    frdiff_db=results.(are).goopgdiff(find(results.(are).goodprdb));
    frdiff_nds=results.(are).goopgdiff(find(results.(are).goodpr));
    frdiff_os=results.(are).goopgdiff(find(results.(are).goodpros));
    frdiff_ob=results.(are).goopgdiff(find(results.(are).goodprob));
    
    frdiff=results.(are).goopgdiff;
    frdiff_recomputed=max_best_plaid_tc-max_best_grating_tc;
    
    % prepare color and select data
    comp_col=[0.2,0.6,0];
    patt_col=[0.9,0.5,0];
    max_best_grating_tc_patt=max_best_grating_tc(pattern_db);
    max_best_grating_tc_comp=max_best_grating_tc(component_db);
    max_best_plaid_tc_patt=max_best_plaid_tc(pattern_db);
    max_best_plaid_tc_comp=max_best_plaid_tc(component_db);
    rel_frd_all=100*(max_best_plaid_tc-max_best_grating_tc)./max_best_tc;
    rel_frd_patt=100*(max_best_plaid_tc(pattern_db)-max_best_grating_tc(pattern_db))./max_best_tc(pattern_db);
    rel_frd_comp=100*(max_best_plaid_tc(component_db)-max_best_grating_tc(component_db))./max_best_tc(component_db);
    frd_all=max_best_plaid_tc-max_best_grating_tc;
    frd_patt=max_best_plaid_tc(pattern_db)-max_best_grating_tc(pattern_db);
    frd_comp=max_best_plaid_tc(component_db)-max_best_grating_tc(component_db);
    
    % compute suppression boolean
    is_suppressed=max_best_grating_tc>max_best_plaid_tc;
    is_suppressed_patt=max_best_grating_tc_patt>max_best_plaid_tc_patt;
    is_suppressed_comp=max_best_grating_tc_comp>max_best_plaid_tc_comp;
    
    n_suppressed=sum(is_suppressed);
    n_suppressed_patt=sum(is_suppressed_patt);
    n_suppressed_comp=sum(is_suppressed_comp);
    
    % observed distribution
    distr_suppressed=[n_suppressed,numel(max_best_grating_tc)-n_suppressed];
    distr_suppressed_patt=[n_suppressed_patt,numel(max_best_grating_tc_patt)-n_suppressed_patt];
    distr_suppressed_comp=[n_suppressed_comp,numel(max_best_grating_tc_comp)-n_suppressed_comp];
    % equivalent uniform distribution
    distr_suppressed_uni_equiv=[sum(distr_suppressed)/2,sum(distr_suppressed)/2];
    distr_suppressed_patt_uni_equiv=[sum(distr_suppressed_patt)/2,sum(distr_suppressed_patt)/2];
    distr_suppressed_comp_uni_equiv=[sum(distr_suppressed_comp)/2,sum(distr_suppressed_comp)/2];
    
    % chi2 different from uniform
    [p_chi2_diff2uni, ~]=chi2homo(distr_suppressed,distr_suppressed_uni_equiv);
    [p_chi2_diff2uni_patt, ~]=chi2homo(distr_suppressed_patt,distr_suppressed_patt_uni_equiv);
    [p_chi2_diff2uni_comp, ~]=chi2homo(distr_suppressed_comp,distr_suppressed_comp_uni_equiv);
    
    % ttest grating plaid response
    [h_tt_all,p_tt_all]=ttest2(max_best_plaid_tc,max_best_grating_tc);
    [h_tt_patt,p_tt_patt]=ttest2(max_best_plaid_tc_patt,max_best_grating_tc_patt);
    [h_tt_comp,p_tt_comp]=ttest2(max_best_plaid_tc_comp,max_best_grating_tc_comp);
    
    % ttest grating plaid response difference
    [h_ttd_all,p_ttd_all]=ttest(frd_all);
    [h_ttd_patt,p_ttd_patt]=ttest(frd_patt);
    [h_ttd_comp,p_ttd_comp]=ttest(frd_comp);
    
    % ttest grating plaid response difference
    [h_ttreld_all,p_ttreld_all]=ttest(rel_frd_all);
    [h_ttreld_patt,p_ttreld_patt]=ttest(rel_frd_patt);
    [h_ttreld_comp,p_ttreld_comp]=ttest(rel_frd_comp);
    
    % get fit coeffs
    [coeffs_all,S_all] = polyfit(max_best_plaid_tc, max_best_grating_tc, 1);
    [coeffs_comp,S_comp] = polyfit(max_best_plaid_tc_comp, max_best_grating_tc_comp, 1);
    [coeffs_patt,S_patt] = polyfit(max_best_plaid_tc_patt, max_best_grating_tc_patt, 1);
    % get fitted values
    fittedX = 0:50;
    fittedY_all = polyval(coeffs_all,fittedX);
    fittedY_comp = polyval(coeffs_comp, fittedX);
    fittedY_patt = polyval(coeffs_patt, fittedX);
    
    % ---------- plot scatter  (responses)   ----------
    
    ff=figure('units','normalized','outerposition',[0 0 1 1]);
    plot(max_best_plaid_tc,max_best_grating_tc,'.','MarkerSize',35,'Color',[0.6,0.6,0.6]); hold on;
    plot(max_best_plaid_tc_patt,max_best_grating_tc_patt,'.','MarkerSize',45,'Color',patt_col)
    plot(max_best_plaid_tc_comp,max_best_grating_tc_comp,'.','MarkerSize',45,'Color',comp_col)
    plot(linspace(0,50,2),...
        linspace(0,50,2),'--','LineWidth',2,'Color',[0,0,0])
    plot(fittedX,fittedY_all,'-','LineWidth',1,'Color',[0.6,0.6,0.6])
    plot(fittedX,fittedY_comp,'-','LineWidth',1,'Color',comp_col)
    plot(fittedX,fittedY_patt,'-','LineWidth',1,'Color',patt_col)
    xlim([0,20]); ylim([0,20]);
    title(['plaid-grating firing rate scatter - ',are],'FontSize',15)
    set(gca,'FontSize',14);
    ylabel('grating peak resonse [Hz]','FontSize',15)
    xlabel('plaid peak resonse [Hz]','FontSize',15)
    ylimit=get(gca,'ylim');
    xlimit=get(gca,'xlim');
    te=text(0.1*xlimit(2),0.95*ylimit(2),['y = ',num2str(coeffs_all(1),'%.2f'),'x + ',num2str(coeffs_all(2),'%.2f')],'FontSize',15,'Color',[0.6,0.6,0.6]);
    te=text(0.1*xlimit(2),0.90*ylimit(2),['y = ',num2str(coeffs_patt(1),'%.2f'),'x + ',num2str(coeffs_patt(2),'%.2f')],'FontSize',15,'Color',patt_col);
    te=text(0.1*xlimit(2),0.85*ylimit(2),['y = ',num2str(coeffs_comp(1),'%.2f'),'x + ',num2str(coeffs_comp(2),'%.2f')],'FontSize',15,'Color',comp_col);
    saveas(ff,['FR_FR_scatter_',are,'.jpg'])
    saveas(ff,['FR_FR_scatter_',are,'.eps'])
    
    % ---------- plot vertical scatter (responses) ----------
    
    ff=figure('units','normalized','outerposition',[0 0 1 1]);
    
    plot(0.5-0.2+0.4*rand(size(max_best_grating_tc)),max_best_grating_tc,'.','MarkerSize',35,'Color',[0.6,0.6,0.6]); hold on;
    plot([0,2],[mean(max_best_grating_tc),mean(max_best_grating_tc)],'LineWidth',2,'LineStyle','--','Color',[0.6,0.6,0.6]); hold on;
    plot(0.5,mean(max_best_grating_tc),'d','MarkerSize',15,'MarkerFaceColor',[0.6,0.6,0.6],'MarkerEdgeColor',[1,1,1]); hold on;
    
    plot(1.5-0.2+0.4*rand(size(max_best_plaid_tc)),max_best_plaid_tc,'.','MarkerSize',35,'Color',[0.6,0.6,0.6]);
    plot([0,2],[mean(max_best_plaid_tc),mean(max_best_plaid_tc)],'LineWidth',2,'LineStyle','--','Color',[0.6,0.6,0.6]); hold on;
    plot(1.5,mean(max_best_plaid_tc),'d','MarkerSize',15,'MarkerFaceColor',[0.6,0.6,0.6],'MarkerEdgeColor',[1,1,1]); hold on;
    
    plot(0.5-0.2+0.4*rand(size(max_best_grating_tc_patt)),max_best_grating_tc_patt,'.','MarkerSize',55,'Color',patt_col); hold on;
    plot([0,2],[mean(max_best_grating_tc_patt),mean(max_best_grating_tc_patt)],'LineWidth',2,'LineStyle','--','Color',patt_col); hold on;
    plot(0.5,mean(max_best_grating_tc_patt),'d','MarkerSize',15,'MarkerFaceColor',patt_col,'MarkerEdgeColor',[1,1,1]); hold on;
    
    plot(0.5-0.2+0.4*rand(size(max_best_grating_tc_comp)),max_best_grating_tc_comp,'.','MarkerSize',55,'Color',comp_col);
    plot([0,2],[mean(max_best_grating_tc_comp),mean(max_best_grating_tc_comp)],'LineWidth',2,'LineStyle','--','Color',comp_col); hold on;
    plot(0.5,mean(max_best_grating_tc_comp),'d','MarkerSize',15,'MarkerFaceColor',comp_col,'MarkerEdgeColor',[1,1,1]); hold on;
    
    plot(1.5-0.2+0.4*rand(size(max_best_plaid_tc_patt)),max_best_plaid_tc_patt,'.','MarkerSize',55,'Color',patt_col); hold on;
    plot([0,2],[mean(max_best_plaid_tc_patt),mean(max_best_plaid_tc_patt)],'LineWidth',2,'LineStyle','--','Color',patt_col); hold on;
    plot(1.5,mean(max_best_plaid_tc_patt),'d','MarkerSize',15,'MarkerFaceColor',patt_col,'MarkerEdgeColor',[1,1,1]); hold on;
    
    plot(1.5-0.2+0.4*rand(size(max_best_plaid_tc_comp)),max_best_plaid_tc_comp,'.','MarkerSize',55,'Color',comp_col);
    plot([0,2],[mean(max_best_plaid_tc_comp),mean(max_best_plaid_tc_comp)],'LineWidth',2,'LineStyle','--','Color',comp_col); hold on;
    plot(1.5,mean(max_best_plaid_tc_comp),'d','MarkerSize',15,'MarkerFaceColor',comp_col,'MarkerEdgeColor',[1,1,1]); hold on;
    
    xlim([0,2])
    ylim([0,20])
    title(['firing rate distributions - ',are],'FontSize',15)
    set(gca,'FontSize',14);
    xticks([0.5,1.5])
    xticklabels({'gratings','plaids'})
    ylabel('peak response [Hz]','FontSize',15)
    xlabel('stimulus identity','FontSize',15)
    ylimit=get(gca,'ylim');
    xlimit=get(gca,'xlim');
    te=text(0.85*xlimit(2),0.95*ylimit(2),['p-tt all = ',num2str(p_tt_all,'%.3f')],'FontSize',15);
    te=text(0.85*xlimit(2),0.90*ylimit(2),['p-tt comp = ',num2str(p_tt_comp,'%.3f')],'FontSize',15);
    te=text(0.85*xlimit(2),0.85*ylimit(2),['p-tt patt = ',num2str(p_tt_patt,'%.3f')],'FontSize',15);
    saveas(ff,['FR_distributions_',are,'.jpg'])
    saveas(ff,['FR_distributions_',are,'.eps'])
    
    % ---------- plot vertical scatter (differences - absolute) ----------
    
    ff=figure('units','normalized','outerposition',[0 0 1 1]);
    
    plot(0.5-0.2+0.4*rand(size(frd_all)),frd_all,'.','MarkerSize',35,'Color',[0.6,0.6,0.6]); hold on;
    plot([0,3],[mean(frd_all),mean(frd_all)],'LineWidth',2,'LineStyle','--','Color',[0.6,0.6,0.6]); hold on;
    plot(0.5,mean(frd_all),'d','MarkerSize',15,'MarkerFaceColor',[0.6,0.6,0.6],'MarkerEdgeColor',[1,1,1]); hold on;
    
    plot(1.5-0.2+0.4*rand(size(frd_comp)),frd_comp,'.','MarkerSize',35,'Color',comp_col);
    plot([0,3],[mean(frd_comp),mean(frd_comp)],'LineWidth',2,'LineStyle','--','Color',comp_col); hold on;
    plot(1.5,mean(frd_comp),'d','MarkerSize',15,'MarkerFaceColor',comp_col,'MarkerEdgeColor',[1,1,1]); hold on;
    
    plot(2.5-0.2+0.4*rand(size(frd_patt)),frd_patt,'.','MarkerSize',55,'Color',patt_col); hold on;
    plot([0,3],[mean(frd_patt),mean(frd_patt)],'LineWidth',2,'LineStyle','--','Color',patt_col); hold on;
    plot(2.5,mean(frd_patt),'d','MarkerSize',15,'MarkerFaceColor',patt_col,'MarkerEdgeColor',[1,1,1]); hold on;
    
    plot([0,3],[0,0],'-','Color',[0,0,0],'LineWidth',2);
    
    xlim([0,3])
    ylim([-10,10])
    title(['firing rate difference distributions - ',are],'FontSize',15)
    set(gca,'FontSize',14);
    xticks([0.5,1.5,2.5])
    xticklabels({'all','component','pattern'})
    ylabel('peak response difference [Hz]','FontSize',15)
    xlabel('cell class','FontSize',15)
    ylimit=get(gca,'ylim');
    xlimit=get(gca,'xlim');
    te=text(0.80*xlimit(2),0.95*ylimit(2),['p-ttd all = ',num2str(p_ttd_all,'%.3f')],'FontSize',15);
    te=text(0.80*xlimit(2),0.90*ylimit(2),['p-ttd comp = ',num2str(p_ttd_comp,'%.3f')],'FontSize',15);
    te=text(0.80*xlimit(2),0.85*ylimit(2),['p-ttd patt = ',num2str(p_ttd_patt,'%.3f')],'FontSize',15);
    te=text(0.80*xlimit(2),0.75*ylimit(2),['mean fr diff all = ',num2str(mean(frd_all),'%.1f'),' Hz'],'FontSize',15);
    te=text(0.80*xlimit(2),0.70*ylimit(2),['mean fr diff comp = ',num2str(mean(frd_comp),'%.1f'),' Hz'],'FontSize',15);
    te=text(0.80*xlimit(2),0.65*ylimit(2),['mean fr diff patt = ',num2str(mean(frd_patt),'%.1f'),' Hz'],'FontSize',15);
    saveas(ff,['FR_difference_distributions_',are,'.jpg'])
    saveas(ff,['FR_difference_distributions_',are,'.eps'])
    
    % ---------- plot vertical scatter (differences - relative) ----------
    
    ff=figure('units','normalized','outerposition',[0 0 1 1]);
    
    plot(0.5-0.2+0.4*rand(size(rel_frd_all)),rel_frd_all,'.','MarkerSize',35,'Color',[0.6,0.6,0.6]); hold on;
    plot([0,3],[mean(rel_frd_all),mean(rel_frd_all)],'LineWidth',2,'LineStyle','--','Color',[0.6,0.6,0.6]); hold on;
    plot(0.5,mean(rel_frd_all),'d','MarkerSize',15,'MarkerFaceColor',[0.6,0.6,0.6],'MarkerEdgeColor',[1,1,1]); hold on;
    
    plot(1.5-0.2+0.4*rand(size(rel_frd_comp)),rel_frd_comp,'.','MarkerSize',35,'Color',comp_col);
    plot([0,3],[mean(rel_frd_comp),mean(rel_frd_comp)],'LineWidth',2,'LineStyle','--','Color',comp_col); hold on;
    plot(1.5,mean(rel_frd_comp),'d','MarkerSize',15,'MarkerFaceColor',comp_col,'MarkerEdgeColor',[1,1,1]); hold on;
    
    plot(2.5-0.2+0.4*rand(size(rel_frd_patt)),rel_frd_patt,'.','MarkerSize',55,'Color',patt_col); hold on;
    plot([0,3],[mean(rel_frd_patt),mean(rel_frd_patt)],'LineWidth',2,'LineStyle','--','Color',patt_col); hold on;
    plot(2.5,mean(rel_frd_patt),'d','MarkerSize',15,'MarkerFaceColor',patt_col,'MarkerEdgeColor',[1,1,1]); hold on;
    
    plot([0,3],[0,0],'-','Color',[0,0,0],'LineWidth',2);
    
    xlim([0,3])
    ylim([-100,100])
    title(['firing rate relative difference distributions - ',are],'FontSize',15)
    set(gca,'FontSize',14);
    xticks([0.5,1.5,2.5])
    xticklabels({'all','component','pattern'})
    ylabel('relative peak response difference (%)','FontSize',15)
    xlabel('cell class','FontSize',15)
    ylimit=get(gca,'ylim');
    xlimit=get(gca,'xlim');
    te=text(0.75*xlimit(2),0.95*ylimit(2),['p-ttreld all = ',num2str(p_ttreld_all,'%.3f')],'FontSize',15);
    te=text(0.75*xlimit(2),0.90*ylimit(2),['p-ttreld comp = ',num2str(p_ttreld_comp,'%.3f')],'FontSize',15);
    te=text(0.75*xlimit(2),0.85*ylimit(2),['p-ttreld patt = ',num2str(p_ttreld_patt,'%.3f')],'FontSize',15);
    te=text(0.75*xlimit(2),0.75*ylimit(2),['mean relative fr diff all = ',num2str(mean(rel_frd_all),'%.1f'),'%'],'FontSize',15);
    te=text(0.75*xlimit(2),0.70*ylimit(2),['mean relative fr diff comp = ',num2str(mean(rel_frd_comp),'%.1f'),'%'],'FontSize',15);
    te=text(0.75*xlimit(2),0.65*ylimit(2),['mean relative fr diff patt = ',num2str(mean(rel_frd_patt),'%.1f'),'%'],'FontSize',15);
    saveas(ff,['FR_relative_difference_distributions_',are,'.jpg'])
    saveas(ff,['FR_relative_difference_distributions_',are,'.eps'])
    
    % ---------- plot suppressed barplot  ----------
    
    ff=figure('units','normalized','outerposition',[0 0 1 1]);
    b=bar([1,2],[100*distr_suppressed/sum(distr_suppressed)],'FaceColor',[0.6,0.6,0.6]); hold on;
    b=bar([4,5],[100*distr_suppressed_comp/sum(distr_suppressed_comp)],'FaceColor',comp_col);
    b=bar([7,8],[100*distr_suppressed_patt/sum(distr_suppressed_patt)],'FaceColor',patt_col);
    xticks([1,2,4,5,7,8]);
    xlim([-1,10])
    ylim([0,110])
    xticklabels({'suppr - all','facil - all','suppr - comp',...
        'facil - comp','suppr - patt','facil - patt'});
    title(['fractions of above diagonal (suppressed) units - ',are])
    ylabel('fraction (%)','FontSize',15)
    set(gca,'FontSize',14);
    ylimit=get(gca,'ylim');
    xlimit=get(gca,'xlim');
    te=text(0.8*xlimit(2),0.95*ylimit(2),['p-chi2uniform all = ',num2str(p_chi2_diff2uni,'%.3f')],'FontSize',15);
    te=text(0.8*xlimit(2),0.90*ylimit(2),['p-chi2uniform comp = ',num2str(p_chi2_diff2uni_comp,'%.3f')],'FontSize',15);
    te=text(0.8*xlimit(2),0.85*ylimit(2),['p-chi2uniform patt = ',num2str(p_chi2_diff2uni_patt,'%.3f')],'FontSize',15);
    saveas(ff,['suppressed_fraction_barplot_',are,'.jpg'])
    saveas(ff,['suppressed_fraction_barplot_',are,'.eps'])
   
end

% cd to old folder
cd(oldfold);
