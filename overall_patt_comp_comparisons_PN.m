% ------------------------- OVERALL PATTERN VERSUS COMPONENT ANALYSIS -------------------------

clear
close all
clc

pars = set_pars_PN();

% set stimulus parameters
listSessions = pars.listSessions;
SF=pars.stimPars.SF;
TF=pars.stimPars.TF;
DIR=pars.stimPars.DIR;
types=pars.neuPars.ctypes;
frame_dur=pars.stimPars.frame_duration;
stim_dur=pars.stimPars.stim_time;
prediction_isi=pars.prediction_isi;

% set analysis parameters
maxlag=pars.STA_depth;
stimwidth=pars.STA_width;
stimheight=pars.STA_height;
interp_factor=pars.interp_factor;
stimlength=pars.stimPars.noiselength;
crop_pixel_size=pars.crop_pixel_size;
contrast_th=pars.contrast_threshold; % 25/07/2017 version used 5.5, now (01/12/2017) 0

% set useful paths
code_folder=pars.code_folder;
addpath(code_folder);
data_folder=pars.processed_data_folder;
addpath(data_folder);
cd(data_folder);
stimulustypes=pars.stimPars.stimulustype;

% load tuning analysis results and indexing
load('Tuning.mat')
load('Indexing.mat')

% create output folder
outfold=['D:\Backups\Personal_bk\PN_acute_analysis\processed_data\rf_quality_comparison_2022_bis'];
outfold=strrep(outfold,'','_');
mkdir(outfold);
oldfold=cd(outfold);

%% perform frame by frame gabor fitting of rfs

% get cell types
celltypes={'component','pattern','unclassified'};
cell_types_codes=[1,2,0];
n_celltypes=numel(celltypes);

% initialize results datastructure
contrast_index_distribs=cell(n_celltypes,3);
selected_cells_idx_distribs=cell(n_celltypes,3);
fitted_rfs_distribs=cell(n_celltypes,3);
fitted_rfs_r2_distribs=cell(n_celltypes,3);

% set options for gabor fitting
gfitoptions.shape = 'elliptical';
gfitoptions.runs  = 5;
gfitoptions.parallel = true;
gfitoptions.visualize = false;

% loop over areas
target_areas={'V1','LM','RL'};

for target_area_idx=1:length(target_areas)
    
    % set area
    are=target_areas{target_area_idx};
    
    % load neuron selection results
    sel_res_file=['D:\Backups\Personal_bk\PN_acute_analysis\processed_data\',are,'_selected_neurons',filesep,'selected_population_',are,'.mat'];
    load(sel_res_file)
    
    % sta analysis results
    sta_res_file=[data_folder,filesep,'RFs_new_datastructure_',are];
    load(sta_res_file)
    
    % loop over cell types
    for cell_type_idx=1:n_celltypes
        
        % select cells of current cell type
        current_cells_idx=find(goomlabel==cell_types_codes(cell_type_idx));
        
        % store current contrast index
        contrast_index_distribs{cell_type_idx,target_area_idx}=gostac(current_cells_idx);
        
        % get current rf
        current_raw_rfs=wsta(:,:,:,current_cells_idx);
        
        % perform rf gabor fitting  for current cell
        current_fitted_rfs=nan(size(current_raw_rfs));
        current_fitted_rfs_r2=nan(size(current_raw_rfs,3),size(current_raw_rfs,4));
        for selected_cell_idx=1:size(current_raw_rfs,4)
            tic
            for frame_idx=1:size(current_raw_rfs,3)
                % fit 2D gabor to current frame
                current_rf_frame=current_raw_rfs(:,:,frame_idx,selected_cell_idx);
                current_rf_frame_resized=imresize(current_rf_frame,2,'bicubic');
                current_rf_frame_resized_fitted = fit2dGabor(current_rf_frame_resized,gfitoptions);
                current_rf_frame_fitted=imresize(current_rf_frame_resized_fitted.patch,1/2,'bicubic');
                current_fitted_rfs(:,:,frame_idx,selected_cell_idx)=current_rf_frame_fitted;
                current_fitted_rfs_r2(frame_idx,selected_cell_idx)=current_rf_frame_resized_fitted.r2;
            end
            toc
            
            % visualize rf gabor fitting results for current cell
            fv=figure('units','normalized','outerposition',[0 0 1 1]);
            v2 = VideoWriter([outfold,filesep,are,'_',celltypes{cell_type_idx},'_n',num2str(selected_cell_idx)],'MPEG-4'); %#ok<TNMLP>            v2.FrameRate = 5;
            v2.FrameRate = 5;
            open(v2);
            for frame_idx=1:size(wsta,3)
                % get current frame to plot
                curr_fr_to_show=size(wsta,3)-frame_idx+1;
                % get max and min of colorscale
                temppixlv=squeeze(current_raw_rfs(:,:,:,selected_cell_idx));
                caxmax=quantile(temppixlv(:),0.99);
                caxmin=quantile(temppixlv(:),0.01);
                % plot sequance of frames
                subplot(1,2,1)
                imagesc(squeeze(current_raw_rfs(:,:,curr_fr_to_show,selected_cell_idx))); colorbar; colormap(gray); caxis([caxmin,caxmax]);
                title('raw RF')
                subplot(1,2,2)
                imagesc(squeeze(current_fitted_rfs(:,:,curr_fr_to_show,selected_cell_idx))); colorbar; colormap(gray); caxis([caxmin,caxmax]);
                title('gabor fitted RF')
                sgtitle([are,' ',celltypes{cell_type_idx},' #',num2str(selected_cell_idx),' - STA frame ',...
                    num2str(curr_fr_to_show),' - r2=  ',num2str(current_fitted_rfs_r2(curr_fr_to_show,selected_cell_idx)),...
                    ' - gosta = ',num2str(contrast_index_distribs{cell_type_idx,target_area_idx}(selected_cell_idx))])
                % append frame to video
                frame = getframe(gcf);
                writeVideo(v2,frame)
            end
            % close video file
            close(v2);
            close all
            
        end
        
        % store current results index
        contrast_index_distribs{cell_type_idx,target_area_idx}=gostac(current_cells_idx);
        fitted_rfs_distribs{cell_type_idx,target_area_idx}=current_fitted_rfs;
        fitted_rfs_r2_distribs{cell_type_idx,target_area_idx}=current_fitted_rfs_r2;
        selected_cells_idx_distribs{cell_type_idx,target_area_idx}=current_cells_idx;
        
    end
    
end

%% perform frame by frame lobe counting of rfs

% initialize results datastructure
contrast_factor_distribs=cell(n_celltypes,3);
relative_area_distribs=cell(n_celltypes,3);
lobe_number_distribs=cell(n_celltypes,3);

% set area limit to consider lobe
input_pars.lobecount_relareath=0.05;
% set binarization thereshold
input_pars.bin_zth=3.5;

% loop over areas
for target_area_idx=1:length(target_areas)
    
    % set area
    are=target_areas{target_area_idx};
    
    % load neuron selection results
    sel_res_file=['D:\Backups\Personal_bk\PN_acute_analysis\processed_data\',are,'_selected_neurons',filesep,'selected_population_',are,'.mat'];
    load(sel_res_file)
    
    % sta analysis results
    sta_res_file=[data_folder,filesep,'RFs_new_datastructure_',are];
    load(sta_res_file)
    
    for cell_type_idx=1:n_celltypes
        
        % select cells of current cell type
        current_cells_idx=find(goomlabel==cell_types_codes(cell_type_idx));
        
        % get current raw rf
        current_raw_Zrfs=Zwsta(:,:,:,current_cells_idx);
        
        % initialize current result vectors
        contrast_factors=nan(size(current_raw_Zrfs,3),size(current_raw_Zrfs,4));
        relative_areas=nan(size(current_raw_Zrfs,3),size(current_raw_Zrfs,4));
        lobe_numbers=nan(size(current_raw_Zrfs,3),size(current_raw_Zrfs,4));
        
        tic
        for selected_cell_idx=1:size(current_raw_Zrfs,4)
            for current_frame_idx=1:size(current_raw_Zrfs,3)
                
                % get and resize current frame
                current_Zrf_frame=squeeze(current_raw_Zrfs(:,:,current_frame_idx,selected_cell_idx));
                current_Zrf_frame_resized=imresize(current_Zrf_frame,10,'bicubic');
                % set cropping function input
                crop_pixel_size=100;
                input_fr=current_Zrf_frame_resized;
                % crop current STA frame
                [ ~, crop_ridx, crop_cidx  ] = apply_crop( input_fr , crop_pixel_size ,0 );
                cropped_frame=input_fr(crop_ridx, crop_cidx);
                % compute STA contrast
                [ ~, ~, ~, contrast_factors(current_frame_idx,selected_cell_idx),...
                    relative_areas(current_frame_idx,selected_cell_idx),...
                    lobe_numbers(current_frame_idx,selected_cell_idx), ~  ] = ...
                    get_shape_params( cropped_frame,input_pars, 0 );
                
            end
        end
        toc
        
        % store current results
        contrast_factor_distribs{cell_type_idx,target_area_idx}=contrast_factors;
        relative_area_distribs{cell_type_idx,target_area_idx}=relative_areas;
        lobe_number_distribs{cell_type_idx,target_area_idx}=lobe_numbers;
        
    end
end

%% recollect OSI and DSI for pattern and component collapsed

% initialize results datastructure
pattern_index_distribs=cell(n_celltypes,3);
Zp_distribs=cell(n_celltypes,3);
Zc_distribs=cell(n_celltypes,3);
osi_distribs=cell(n_celltypes,3);
dsi_distribs=cell(n_celltypes,3);

% loop over areas
for target_area_idx=1:length(target_areas)
    
    % set area
    are=target_areas{target_area_idx};
    
    % load neuron selection results
    sel_res_file=[are,'_selected_neurons',filesep,'selected_population_',are,'.mat'];
    load(sel_res_file)
    
    %     % recompute pattern and componant cell classification
    %     bool_above_zp=(goozp-max(goozc,zeros(size(goozc))))>=1.28;
    %     bool_above_zc=(goozc-max(goozp,zeros(size(goozp))))>=1.28;
    %     ctindex=zeros(size(bool_above_zp));
    %     ctindex(bool_above_zc)=1;
    %     ctindex(bool_above_zp)=2;
    
    for cell_type_idx=1:n_celltypes
        tic
        
        % select cells of current cell type
        current_cells_idx=find(goomlabel==cell_types_codes(cell_type_idx));
        
        % get current neurons OSI value
        osi_distribs{cell_type_idx,target_area_idx}=gooosi(current_cells_idx);
        % get current neurons DSI value
        dsi_distribs{cell_type_idx,target_area_idx}=goodsi(current_cells_idx);
        % get current neurons PI value
        pattern_index_distribs{cell_type_idx,target_area_idx}=goopi(current_cells_idx);
        % get current neurons Zp value
        Zp_distribs{cell_type_idx,target_area_idx}=goozp(current_cells_idx);
        % get current neurons Zc value
        Zc_distribs{cell_type_idx,target_area_idx}=goozc(current_cells_idx);
        
        toc
    end
end

%% collapse areas together for plotting

% set threshold contsrast factor to select  area
cf_th_to_select_area=10;

% initialize collapsed variables
collapsed_log_contrast_index_distribs=cell(n_celltypes,1);
collapsed_max_fitted_rfs_r2_distribs=cell(n_celltypes,1);
collapsed_unrl_fitted_rfs_r2_distribs=cell(n_celltypes,1);
collapsed_unrl_contrast_factor_distribs=cell(n_celltypes,1);
collapsed_unrl_relative_area_distribs=cell(n_celltypes,1);
collapsed_unrl_relative_area_distribs_unfilt=cell(n_celltypes,1);
collapsed_unrl_lobe_number_distribs=cell(n_celltypes,1);
collapsed_best_lobe_number_distribs=cell(n_celltypes,1);
collapsed_pattern_index_distribs=cell(n_celltypes,1);
collapsed_Zp_distribs=cell(n_celltypes,1);
collapsed_Zc_distribs=cell(n_celltypes,1);
collapsed_osi_distribs=cell(n_celltypes,1);
collapsed_dsi_distribs=cell(n_celltypes,1);
% loop over cell types
for cell_type_idx=1:n_celltypes
    % loop over areas
    for target_area_idx=1:length(target_areas)
        tic
        
        % rf gabor fit results
        collapsed_log_contrast_index_distribs{cell_type_idx}=cat(1,...
            collapsed_log_contrast_index_distribs{cell_type_idx},...
            log10(contrast_index_distribs{cell_type_idx,target_area_idx}));
        collapsed_max_fitted_rfs_r2_distribs{cell_type_idx}=cat(1,...
            collapsed_max_fitted_rfs_r2_distribs{cell_type_idx},...
            max(fitted_rfs_r2_distribs{cell_type_idx,target_area_idx})');
        collapsed_unrl_fitted_rfs_r2_distribs{cell_type_idx}=cat(1,...
            collapsed_unrl_fitted_rfs_r2_distribs{cell_type_idx},...
            (fitted_rfs_r2_distribs{cell_type_idx,target_area_idx}(:)));
        % rf shape analysis results
        collapsed_unrl_contrast_factor_distribs{cell_type_idx}=cat(1,...
            collapsed_unrl_contrast_factor_distribs{cell_type_idx},...
            (contrast_factor_distribs{cell_type_idx,target_area_idx}(:)));
        temp_cf=contrast_factor_distribs{cell_type_idx,target_area_idx}(:);
        temp_a=relative_area_distribs{cell_type_idx,target_area_idx}(:);
        collapsed_unrl_relative_area_distribs{cell_type_idx}=cat(1,...
            collapsed_unrl_relative_area_distribs{cell_type_idx},...
            (temp_a(temp_cf>=cf_th_to_select_area)));
        collapsed_unrl_relative_area_distribs_unfilt{cell_type_idx}=cat(1,...
            collapsed_unrl_relative_area_distribs_unfilt{cell_type_idx},...
            (temp_a));
        collapsed_unrl_lobe_number_distribs{cell_type_idx}=cat(1,...
            collapsed_unrl_lobe_number_distribs{cell_type_idx},...
            (lobe_number_distribs{cell_type_idx,target_area_idx}(:)));
        [~,best_idx]=max(contrast_factor_distribs{cell_type_idx,target_area_idx},[],1);
        tempvec=nan(1,numel(best_idx));
        for kk=1:numel(best_idx)
            tempvec(kk)=lobe_number_distribs{cell_type_idx,target_area_idx}(best_idx(kk),kk);
        end
        collapsed_best_lobe_number_distribs{cell_type_idx}=cat(1,...
            collapsed_best_lobe_number_distribs{cell_type_idx},...
            tempvec');
        % selectivity indexes
        collapsed_pattern_index_distribs{cell_type_idx}=cat(1,...
            collapsed_pattern_index_distribs{cell_type_idx},...
            (pattern_index_distribs{cell_type_idx,target_area_idx}));
        collapsed_Zp_distribs{cell_type_idx}=cat(1,...
            collapsed_Zp_distribs{cell_type_idx},...
            (Zp_distribs{cell_type_idx,target_area_idx}));
        collapsed_Zc_distribs{cell_type_idx}=cat(1,...
            collapsed_Zc_distribs{cell_type_idx},...
            (Zc_distribs{cell_type_idx,target_area_idx}));
        temposi=osi_distribs{cell_type_idx,target_area_idx};
        temposi(temposi<0)=NaN;
        collapsed_osi_distribs{cell_type_idx}=cat(1,...
            collapsed_osi_distribs{cell_type_idx},...
            (temposi));
        tempdsi=dsi_distribs{cell_type_idx,target_area_idx};
        tempdsi(tempdsi<0)=NaN;
        collapsed_dsi_distribs{cell_type_idx}=cat(1,...
            collapsed_dsi_distribs{cell_type_idx},...
            (tempdsi));
        
        toc
    end
end

%% save rf fits and lobe counting results

% save results
save([outfold,filesep,'rf_fitting_results.mat'],...
    'contrast_index_distribs',...
    'fitted_rfs_distribs',...
    'fitted_rfs_r2_distribs',...
    'selected_cells_idx_distribs',...
    'pattern_index_distribs',...
    'Zp_distribs',...
    'Zc_distribs',...
    'osi_distribs',...
    'dsi_distribs',...
    'collapsed_log_contrast_index_distribs',...
    'collapsed_max_fitted_rfs_r2_distribs',...
    'collapsed_unrl_fitted_rfs_r2_distribs',...
    'collapsed_unrl_contrast_factor_distribs',...
    'collapsed_unrl_relative_area_distribs',...
    'collapsed_unrl_lobe_number_distribs',...
    'collapsed_best_lobe_number_distribs',...
    'collapsed_pattern_index_distribs',...
    'collapsed_Zp_distribs',...
    'collapsed_Zc_distribs',...
    'collapsed_osi_distribs',...
    'collapsed_dsi_distribs'...
    )

load([outfold,filesep,'rf_fitting_results.mat'])

%% plot gabor fit goodness of fit distributions violins

% decide wheter to use max or unrolled
distribtouse_orig=collapsed_unrl_fitted_rfs_r2_distribs([1,2]); % collapsed_max_fitted_rfs_r2_distribs
distribtouse=cell(size(distribtouse_orig));
% filter according to selectivity
th_val=-1;
for ii=1:numel(distribtouse_orig)
    temp=repmat(collapsed_dsi_distribs{ii},[1,10])';
    distribtouse_filtvals=temp(:);
    distribtouse{ii}=distribtouse_orig{ii}(distribtouse_filtvals>=th_val);
end

% initialize figure
fighand1=figure('units','normalized','outerposition',[0 0 1 1]);
inputpars.inputaxh=gca;
hold(inputpars.inputaxh,'on')
% set settings for violin distribution plotting
inputpars.boxplotwidth=0.4;%0.5;
inputpars.boxplotlinewidth=2;
inputpars.densityplotwidth=0.4;%0.5;
inputpars.yimtouse=[-0.1,0.9];
% inputpars.yimtouse=[0,8];
inputpars.scatterjitter=inputpars.boxplotlinewidth*0.1;
inputpars.scatteralpha=0.25;
inputpars.scattersize=20;
inputpars.distralpha=0.5;
inputpars.xlabelstring=[];
inputpars.ylabelstring='goodness of Gabor fit (r2)';
inputpars.titlestring=['goodness of Gabor fit - all areas ( #component fr = ',...
    num2str(numel(distribtouse{1}),'%.0f'),...
    ' - #pattern fr = ',num2str(numel(distribtouse{2}),'%.0f'),' )'];
inputpars.boolscatteron=1;
inputpars.ks_bandwidth=0.045; % inputpars.ks_bandwidth=0.25;
inputpars.xlimtouse=[-0.5,3.5]; %[-1,5];
% plot violins
inputadata.inputdistrs=distribtouse;
inputpars.n_distribs=numel(inputadata.inputdistrs);
inputpars.dirstrcenters=(1:inputpars.n_distribs);
inputpars.xtickslabelvector={'component','pattern'};
inputpars.distrcolors{1}=[50,200,0]./255;
inputpars.distrcolors{2}=[255,150,0]./255;
inputaxh = plot_violinplot_PN_new(inputadata,inputpars);
pvalw = ranksum(distribtouse{1},distribtouse{2});
[junk,pvalt] = ttest2(distribtouse{1},distribtouse{2});
text(-0.4,0.5,['median diff p = ',num2str(pvalw)],'fontsize',12);
text(-0.4,0.45,['mean diff p = ',num2str(pvalt)],'fontsize',12);
set(gca,'fontsize',12)
saveas(fighand1,[outfold,filesep,'patt_comp_gabor_gof_all_areas'],'jpg')
print(fighand1,'-depsc','-painters',[[outfold,filesep,'patt_comp_gabor_gof_all_areas'],'.eps'])

%% plot shape features distributions violins

distrtoplotlist_orig{1}=collapsed_unrl_contrast_factor_distribs;
distrtoplotlist_orig{2}=collapsed_unrl_relative_area_distribs;
distrtoplotlist_orig{3}=collapsed_best_lobe_number_distribs;
distrtoplotlist{1}=cell(size(collapsed_unrl_contrast_factor_distribs));
distrtoplotlist{2}=distrtoplotlist_orig{2};
distrtoplotlist{3}=cell(size(collapsed_best_lobe_number_distribs));
% filter according to selectivity
th_val=0.33;
for ii=1:numel(distrtoplotlist_orig{1})
    temp=repmat(collapsed_dsi_distribs{ii},[1,10])';
    distribtouse_filtvals=temp(:);
    distrtoplotlist{1}{ii}=distrtoplotlist_orig{1}{ii}(distribtouse_filtvals>=th_val);
end
for ii=1:numel(distrtoplotlist_orig{3})
    temp=repmat(collapsed_dsi_distribs{ii},[1,1])';
    distribtouse_filtvals=temp(:);
    distrtoplotlist{3}{ii}=distrtoplotlist_orig{3}{ii}(distribtouse_filtvals>=th_val);
end

ylabellist{1}='contrast factor';
ylabellist{2}='relative area';
yimtouselist{1}=[2.5,20];
yimtouselist{2}=[0.025,0.35];
ks_ban{1}=0.5;
ks_ban{2}=0.01;
% initialize figure
fighand2=figure('units','normalized','outerposition',[0 0 1 1]);
for jj=1:2
    % decide wheter to use max or unrolled
    distribtouse=distrtoplotlist{jj}; % collapsed_max_fitted_rfs_r2_distribs
    subplot(1,3,jj)
    inputpars.inputaxh=gca;
    hold(inputpars.inputaxh,'on')
    % set settings for violin distribution plotting
    inputpars.boxplotwidth=0.4;%0.5;
    inputpars.boxplotlinewidth=2;
    inputpars.densityplotwidth=0.4;%0.5;
    inputpars.yimtouse=yimtouselist{jj};
    % inputpars.yimtouse=[0,8];
    inputpars.scatterjitter=inputpars.boxplotlinewidth*0.1;
    inputpars.scatteralpha=0.15;
    inputpars.scattersize=20;
    inputpars.distralpha=0.5;
    inputpars.xlabelstring=[];
    inputpars.ylabelstring=ylabellist{jj};
    inputpars.titlestring=[ylabellist{jj},' ( #comp fr = ',...
        num2str(numel(distribtouse{1}),'%.0f'),...
        ' - #patt fr = ',num2str(numel(distribtouse{2}),'%.0f'),' )'];
    inputpars.boolscatteron=1;
    inputpars.ks_bandwidth=ks_ban{jj}; % inputpars.ks_bandwidth=0.25;
    inputpars.xlimtouse=[-0.5,4.5]; %[-1,5];
    % plot violins
    inputadata.inputdistrs=distribtouse;
    inputpars.n_distribs=numel(inputadata.inputdistrs);
    inputpars.dirstrcenters=(1:inputpars.n_distribs);
    inputpars.xtickslabelvector={'component','pattern','unclassified'};
    inputpars.distrcolors{1}=[50,200,0]./255;
    inputpars.distrcolors{2}=[255,150,0]./255;
    inputpars.distrcolors{3}=[0,0,0]./(3*255);
    inputaxh = plot_violinplot_PN_new(inputadata,inputpars);
    pvalw = ranksum(distribtouse{1},distribtouse{2});
    [junk,pvalt] = ttest2(distribtouse{1},distribtouse{2});
    text(-0.2,0.95*max(yimtouselist{jj}),['median diff p = ',num2str(pvalw)],'fontsize',12);
    text(-0.2,0.92*max(yimtouselist{jj}),['mean diff p = ',num2str(pvalt)],'fontsize',12);
    set(gca,'fontsize',12)
end
subplot(1,3,3)
comp_lobe_n_dis=distrtoplotlist{3}{1};
patt_lobe_n_dis=distrtoplotlist{3}{2};
temp=[comp_lobe_n_dis;patt_lobe_n_dis];
xtouse=[min(temp):1:max(temp)];
[N_comp_lobe_n,~] = hist(comp_lobe_n_dis,xtouse); %#ok<HIST>
p_comp_lobe_n=N_comp_lobe_n./sum(N_comp_lobe_n);
[N_patt_lobe_n,X] = hist(patt_lobe_n_dis,xtouse); %#ok<HIST>
p_patt_lobe_n=N_patt_lobe_n./sum(N_patt_lobe_n);
hold on;
bar(10*X,p_comp_lobe_n,...
    'facecolor',inputpars.distrcolors{1},...
    'edgecolor',[0,0,0],...
    'facealpha',0.5,...
    'linewidth',1.5,...
    'barwidth',0.3)
bar(10*X+3,p_patt_lobe_n,...
    'facecolor',inputpars.distrcolors{2},...
    'edgecolor',[0,0,0],...
    'linewidth',1.5,...
    'facealpha',0.5,...
    'barwidth',0.3)
titlestring=['lobe count',' ( #comp fr = ',...
    num2str(numel(distrtoplotlist{3}{1}),'%.0f'),...
    ' - #patt fr = ',num2str(numel(distrtoplotlist{3}{2}),'%.0f'),' )'];
title(titlestring)
ylim([0,1])
[chi2stat,p_chitest] = chiSquareTest([N_comp_lobe_n;N_patt_lobe_n]);
text(20,0.92,['chi square p = ',num2str(p_chitest)],'fontsize',12);
xticks(10*X+1.5)
xticklabels(X)
xlabel('# of lobes')
ylabel('fraction of cells')
set(gca,'fontsize',12)
sgtitle(['RF shape analysis - all areas'])
saveas(fighand2,[outfold,filesep,'patt_comp_rf_shape_analysis_all_areas_DS'],'jpg')
print(fighand2,'-depsc','-painters',[[outfold,filesep,'patt_comp_rf_shape_analysis_all_areas_DS'],'.eps'])

%% plot selectivity index violins

distrtoplotlist{1}=collapsed_osi_distribs;
distrtoplotlist{2}=collapsed_dsi_distribs;
distrtoplotlist{3}=collapsed_pattern_index_distribs;
ylabellist{1}='OSI';
ylabellist{2}='DSI';
ylabellist{3}='PI';
yimtouselist{1}=[-0.2,1.2];
yimtouselist{2}=[-0.2,1.2];
yimtouselist{3}=[-7,7];
ks_ban{1}=0.07;
ks_ban{2}=0.07;
ks_ban{3}=0.3;
% initialize figure
fighand3=figure('units','normalized','outerposition',[0 0 1 1]);
for jj=1:n_celltypes
    % decide wheter to use max or unrolled
    distribtouse=distrtoplotlist{jj}; % collapsed_max_fitted_rfs_r2_distribs
    subplot(1,3,jj)
    inputpars.inputaxh=gca;
    hold(inputpars.inputaxh,'on')
    % set settings for violin distribution plotting
    inputpars.boxplotwidth=0.4;%0.5;
    inputpars.boxplotlinewidth=2;
    inputpars.densityplotwidth=0.4;%0.5;
    inputpars.yimtouse=yimtouselist{jj};
    % inputpars.yimtouse=[0,8];
    inputpars.scatterjitter=inputpars.boxplotlinewidth*0.1;
    inputpars.scatteralpha=0.15;
    inputpars.scattersize=20;
    inputpars.distralpha=0.5;
    inputpars.xlabelstring=[];
    inputpars.ylabelstring=ylabellist{jj};
    inputpars.titlestring=[ylabellist{jj},' ( #comp fr = ',...
        num2str(numel(distribtouse{1}),'%.0f'),...
        ' - #patt fr = ',num2str(numel(distribtouse{2}),'%.0f'),' )'];
    inputpars.boolscatteron=1;
    inputpars.ks_bandwidth=ks_ban{jj}; % inputpars.ks_bandwidth=0.25;
    inputpars.xlimtouse=[-0.5,4.5]; %[-1,5];
    % plot violins
    inputadata.inputdistrs=distribtouse;
    inputpars.n_distribs=numel(inputadata.inputdistrs);
    inputpars.dirstrcenters=(1:inputpars.n_distribs);
    inputpars.xtickslabelvector={'component','pattern'};
    inputpars.distrcolors{1}=[50,200,0]./255;
    inputpars.distrcolors{2}=[255,150,0]./255;
    inputaxh = plot_violinplot_PN_new(inputadata,inputpars);
    pvalw = ranksum(distribtouse{1},distribtouse{2});
    [junk,pvalt] = ttest2(distribtouse{1},distribtouse{2});
    if jj==3
        hold on;
        plot([min(inputpars.xlimtouse), max(inputpars.xlimtouse)],[0,0],'--','linewidth',2,'color',[0.5,0.5,0.5])
    else
        text(-0.2,0.95*max(yimtouselist{jj}),['median diff p = ',num2str(pvalw)],'fontsize',12);
        text(-0.2,0.92*max(yimtouselist{jj}),['mean diff p = ',num2str(pvalt)],'fontsize',12);
    end
    set(gca,'fontsize',12)
end
sgtitle(['Selectivity analysis - all areas'])
saveas(fighand3,[outfold,filesep,'patt_comp_selectivity_analysis_all_areas'],'jpg')

%% plot partial correlation plot - all areas

fighand4=figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
npatt_2=0;
ncomp_2=0;
for cell_class_idx=1:n_celltypes
    scatter(collapsed_Zc_distribs{cell_class_idx},collapsed_Zp_distribs{cell_class_idx},100,...
        'MarkerFaceColor',inputpars.distrcolors{cell_class_idx},'MarkerEdgeColor',inputpars.distrcolors{cell_class_idx},...
        'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
    valid_neu_idx=collapsed_dsi_distribs{cell_class_idx}>=0.33;
    tempZp=collapsed_Zp_distribs{cell_class_idx};
    tempZc=collapsed_Zc_distribs{cell_class_idx};
    plot(tempZc(valid_neu_idx),tempZp(valid_neu_idx),'.','MarkerSize',25,'Color',inputpars.distrcolors{cell_class_idx});
    if cell_class_idx==1
        ncomp_2=sum(valid_neu_idx);
    else
        npatt_2=sum(valid_neu_idx);
    end
end
npatt_1=numel(collapsed_Zc_distribs{2});
ncomp_1=numel(collapsed_Zc_distribs{1});
line([0 6], [1.28 7.28],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
line([1.28 7.28], [0 6],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
line([1.28 1.28], [-4 0],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
line([-4 0], [1.28 1.28],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
xlabel('Zc'); ylabel('Zp');
title(['Partial correlation scatter - all areas']);
text(-2.2,0.95*6.5,['n patt (DSI>0) = ',num2str(npatt_1)],'fontsize',12);
text(-2.2,0.92*6.5,['n comp (DSI>0) = ',num2str(ncomp_1)],'fontsize',12);
text(-2.2,0.89*6.5,['n patt (DSI>0.33) = ',num2str(npatt_2)],'fontsize',12);
text(-2.2,0.86*6.5,['n comp (DSI>0.33) = ',num2str(ncomp_2)],'fontsize',12);
set(gca,'fontsize',12);
ylim(gca,[-3.5,6.5])
xlim(gca,[-3.5,6.5])
axis square
saveas(fighand4,[outfold,filesep,'patt_comp_partial_correlation_all_areas'],'jpg')
print(fighand4,'-depsc','-painters',[[outfold,filesep,'patt_comp_partial_correlation_all_areas'],'.eps'])

%% plot partial correlation plot - each area separtely

areas={'V1','LM','RL'};
fighand4=figure('units','normalized','outerposition',[0 0 1 1]);
for area_idx=1:3
    subplot(1,3,area_idx)
    hold on;
    npatt_2=0;
    ncomp_2=0;
    for cell_class_idx=1:n_celltypes
        scatter(Zc_distribs{cell_class_idx,area_idx},Zp_distribs{cell_class_idx,area_idx},100,...
            'MarkerFaceColor',inputpars.distrcolors{cell_class_idx},'MarkerEdgeColor',inputpars.distrcolors{cell_class_idx},...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
        valid_neu_idx=dsi_distribs{cell_class_idx,area_idx}>=0.33;
        tempZp=Zp_distribs{cell_class_idx,area_idx};
        tempZc=Zc_distribs{cell_class_idx,area_idx};
        plot(tempZc(valid_neu_idx),tempZp(valid_neu_idx),'.','MarkerSize',25,'Color',inputpars.distrcolors{cell_class_idx});
        if cell_class_idx==1
            ncomp_2=sum(valid_neu_idx);
        else
            npatt_2=sum(valid_neu_idx);
        end
    end
    npatt_1=numel(Zc_distribs{2,area_idx});
    ncomp_1=numel(Zc_distribs{1,area_idx});
    line([0 6], [1.28 7.28],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
    line([1.28 7.28], [0 6],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
    line([1.28 1.28], [-4 0],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
    line([-4 0], [1.28 1.28],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
    xlabel('Zc'); ylabel('Zp');
    title(['Partial correlation scatter - ',areas{area_idx}]);
    text(-2.2,0.95*6.5,['n patt (DSI>0) = ',num2str(npatt_1)],'fontsize',12);
    text(-2.2,0.90*6.5,['n comp (DSI>0) = ',num2str(ncomp_1)],'fontsize',12);
    text(-2.2,0.85*6.5,['n patt (DSI>0.33) = ',num2str(npatt_2)],'fontsize',12);
    text(-2.2,0.80*6.5,['n comp (DSI>0.33) = ',num2str(ncomp_2)],'fontsize',12);
    set(gca,'fontsize',12);
    ylim(gca,[-3.5,6.5])
    xlim(gca,[-3.5,6.5])
    axis square
end
% saveas(fighand4,[outfold,filesep,'patt_comp_partial_correlation_all_areas_separately'],'jpg')
% print(fighand4,'-depsc','-painters',[[outfold,filesep,'patt_comp_partial_correlation_all_areas_separately'],'.eps'])