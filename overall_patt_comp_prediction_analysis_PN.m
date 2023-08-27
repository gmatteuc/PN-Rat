% ------------------------- OVERALL PATTERN VERSUS COMPONENT PREDICTION ANALYSIS -------------------------

clear
close all
clc

% set pars
pars = set_pars_PN();

% add home code path
addpath('D:\Backups\Personal_bk\PN_acute_analysis\scripts\new_analyses_2022')

% set cett type label to use
labelstouse='oldbest'; %'old' 'new' 'oldbest' 'oldbest_intersection'

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
% crop_pixel_size=pars.crop_pixel_size;
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
outfold=['D:\Backups\Personal_bk\PN_acute_analysis\processed_data\rf_quality_comparison_2022_final_recovery_test'];
outfold=strrep(outfold,'','_');
mkdir(outfold);
oldfold=cd(outfold);
addpath(outfold)

%% perform prediction anlysis with fitted sta filters

% set r2 threshold for prediction analysis
r2_th_for_pred_anlaysis=0;%0.3;

% set video framerate
videoframerate=30;

% get contrast index by cell type and area
cell_types={'component','pattern','unclassified'};

% loop over areas
target_areas={'V1','LM','RL'};

% if needed rerun prediction
bool_redo_predictions=0;
if bool_redo_predictions
    
    % set neuron idx to plot NB: added for plotting
    % neuron_idx_toplot={[44],[403,396,438],[625]};
    % neuron_idx_toplot={[44,137,890,532,547],[427,403,396,438],[625]};
    
    % loop over areas
    for target_area_idx=1:length(target_areas)%1:length(target_areas)
        
        % set area
        are=target_areas{target_area_idx};
        
        % get current neuron idx to plot NB: added for plotting
        %     current_neuron_idx_toplot=neuron_idx_toplot{target_area_idx};
        
        % loop over cell types
        for cell_type_idx=1:length(cell_types)
            
            % set cell type
            cty=cell_types{cell_type_idx};
            
            % load neuron selection results
            sel_res_file=['D:\Backups\Personal_bk\PN_acute_analysis\processed_data',filesep,are,'_selected_neurons',filesep,'selected_population_',are,'.mat'];
            load(sel_res_file)
            
            % load fitted sta analysis results
            sta_res_file=[outfold,filesep,'rf_fitting_results.mat'];
            load(sta_res_file)
            
            % select current category neurons with good enough Gabor fitted rf
            current_r2s_per_frame=fitted_rfs_r2_distribs{cell_type_idx,target_area_idx};
            valid_neurons_idx=find(max(current_r2s_per_frame,[],1)>r2_th_for_pred_anlaysis);
            current_selected_neu_idx=selected_cells_idx_distribs{cell_type_idx,target_area_idx}(valid_neurons_idx);
            selectedsi_to_use=selectedsi(current_selected_neu_idx); % inex within the ogiginal good .. structure
            selectedci_to_use=contrast_index_distribs{cell_type_idx,target_area_idx}(valid_neurons_idx);
            %filters_to_use=fitted_rfs_distribs{cell_type_idx,target_area_idx}(:,:,:,valid_neurons_idx);
            r2s_per_frame_to_use=current_r2s_per_frame(:,valid_neurons_idx);
            
            % recover original STA for display
            sta_res_file=[data_folder,filesep,'RFs_new_datastructure_',are];
            load(sta_res_file)
            cell_types_codes=[1,2,0];
            current_cells_idx=find(goomlabel==cell_types_codes(cell_type_idx));
            filters_to_use_raw=Zwsta(:,:,:,current_cells_idx(valid_neurons_idx));
            cis_per_frame_to_use=contrast(:,current_cells_idx(valid_neurons_idx));
            filters_to_use=Zwsta(:,:,:,current_cells_idx(valid_neurons_idx));%fitted_rfs_distribs{cell_type_idx,target_area_idx}(:,:,:,valid_neurons_idx);
            
            % initialize prediction storage variables
            stimulus_binnumber=floor(stim_dur/frame_dur);
            isi_binnumber=prediction_isi;
            prediction_binnumber=isi_binnumber + stimulus_binnumber + isi_binnumber;
            pred_FR=zeros(videoframerate,length(DIR),length(selectedsi_to_use),numel(stimulustypes));
            obs_FR=zeros(videoframerate,length(DIR),length(selectedsi_to_use),numel(stimulustypes));
            pred_COUNT=zeros(length(DIR),length(selectedsi_to_use),numel(stimulustypes));
            obs_COUNT=zeros(length(DIR),length(selectedsi_to_use),numel(stimulustypes));
            explvar_sD=zeros(length(selectedsi_to_use),numel(stimulustypes));
            explvar_D=zeros(length(selectedsi_to_use),numel(stimulustypes));
            explvar_T=zeros(length(selectedsi_to_use),numel(stimulustypes));
            pDIRdiff=zeros(length(selectedsi_to_use),numel(stimulustypes));
            pORIdiff=zeros(length(selectedsi_to_use),numel(stimulustypes));
            PI_O=zeros(length(selectedsi_to_use),1);
            Zp_O=zeros(length(selectedsi_to_use),1);
            Zc_O=zeros(length(selectedsi_to_use),1);
            Rp_O=zeros(length(selectedsi_to_use),1);
            Rc_O=zeros(length(selectedsi_to_use),1);
            PI_P=zeros(length(selectedsi_to_use),1);
            Zp_P=zeros(length(selectedsi_to_use),1);
            Zc_P=zeros(length(selectedsi_to_use),1);
            Rp_P=zeros(length(selectedsi_to_use),1);
            Rc_P=zeros(length(selectedsi_to_use),1);
            ctlabel_O=zeros(length(selectedsi_to_use),1);
            ctlabel_P=zeros(length(selectedsi_to_use),1);
            obs_OSI=zeros(length(selectedsi_to_use),1);
            obs_DSI=zeros(length(selectedsi_to_use),1);
            alphaopt=zeros(length(selectedsi_to_use),1);
            betaopt=zeros(length(selectedsi_to_use),1);
            good_contrast=zeros(length(selectedsi_to_use),1);
            good_power=zeros(length(selectedsi_to_use),1);
            
            % loop over selected neurons
            for i=1:length(selectedsi_to_use)
                %            if sum(current_neuron_idx_toplot==selectedsi_to_use{i}(1)) % NB: added for plotting
                
                % get current neuron ------------------------------------------------------------------
                
                % get back original good neuron index os and dsi
                orig_i=current_selected_neu_idx(i);
                obs_DSI(i)=goodsi(orig_i);
                obs_OSI(i)=goodsi(orig_i);
                
                % get back STA filter
                selectedSTA=filters_to_use(:,:,:,i);
                selectedSTA_raw=filters_to_use_raw(:,:,:,i); % NB: added for plotting
                n=selectedsi_to_use{i}(1);
                stimtoplot=cell(1,numel(stimulustypes));
                rastertoplot=cell(1,numel(stimulustypes));
                pdirtoplot=cell(1,numel(stimulustypes));
                psftoplot=cell(1,numel(stimulustypes));
                ptftoplot=cell(1,numel(stimulustypes));
                
                % loop over stimulus types
                for k=1:numel(stimulustypes)
                    
                    % get preferred TF condition
                    pSF=SF(1); %SF(1);
                    pTF=TF(selectedsi_to_use{i}(3));
                    [~,idxmax]=randmax(tuning_curve(:,n,pSF==SF,pTF==TF,k));
                    pDIR=DIR(idxmax);
                    % get non-preferred TF condition
                    npSF=SF(1); %SF(1);
                    npTF=TF(setdiff(1:2,selectedsi_to_use{i}(3)));
                    [~,idxmax]=randmax(tuning_curve(:,n,pSF==SF,pTF==TF,k));
                    npDIR=DIR(idxmax);
                    
                    % collect p dir to plot
                    pdirtoplot{k}=pDIR;
                    psftoplot{k}=pSF;
                    ptftoplot{k}=pTF;
                    
                    % prepare inputs for nonlinearity optimization (training) ------------------------------------------------------------------
                    if k==1
                        tic
                        
                        % set parameters for predictions
                        input_interp_factor=interp_factor;
                        input_sf=npSF;
                        input_tf=npTF;
                        input_dirs=DIR;
                        input_stimtype=stimulustypes{k};
                        input_r2_per_frame=r2s_per_frame_to_use(:,i);
                        input_STA=filters_to_use(:,:,:,i);
                        input_neuron_num=n;
                        %NB: shorten this time loading stimuli all at a time
                        [input_PSTHs, ~, ~, input_F, input_S] = get_inputs_for_neuron_responses_prediction(...
                            input_interp_factor,input_neuron_num,input_sf,input_tf,input_dirs,input_stimtype,input_r2_per_frame,input_STA,r2_th_for_pred_anlaysis);
                        % output progress message
                        fprintf(['collecting inputs for ',input_stimtype,...
                            's ( SF=',num2str(pSF),' TF=',num2str(pTF),' ) neuron ',num2str(n),' collected ...\n']);
                        
                        toc
                        
                        % optimize nonlinearity (training) ------------------------------------------------------------------
                        
                        tic
                        
                        % set initial conditions
                        sr=videoframerate;
                        countingwindowlim=[0,1];
                        alpha0=1;
                        beta0=1;
                        % optimize nonlinearity parameter
                        costfunc = @(beta) prediction_error(input_PSTHs,input_S,input_dirs,...
                            pDIR,alpha0,beta,input_F,sr,countingwindowlim);
                        cost0=costfunc(beta0);
                        lb = [1]; ub = [2]; A = []; b = []; Aeq = []; beq = []; %#ok<NBRAK>
                        [temp_betaopt,fval,exitflag,output] = fmincon(costfunc,beta0,A,b,Aeq,beq,lb,ub);
                        %                     [temp_betaopt,fval,exitflag,output]=fminunc(costfunc,beta0);
                        alphaopt(i)=alpha0;
                        betaopt(i)=temp_betaopt;
                        toc
                        
                    end
                    
                    % prepare inputs for prediction of tuning curve and psth (testing) ------------------------------------------------------------------
                    
                    tic
                    
                    % set parameters for predictions
                    input_interp_factor=interp_factor;
                    input_sf=pSF;
                    input_tf=pTF;
                    input_dirs=DIR;
                    input_stimtype=stimulustypes{k};
                    input_r2_per_frame=r2s_per_frame_to_use(:,i);
                    input_STA=filters_to_use(:,:,:,i);
                    input_neuron_num=n;
                    [input_PSTHs, input_RASTERs, input_STIMULIs, input_F, input_S] = get_inputs_for_neuron_responses_prediction(...
                        input_interp_factor,input_neuron_num,input_sf,input_tf,input_dirs,input_stimtype,input_r2_per_frame,input_STA,r2_th_for_pred_anlaysis);
                    % output progress message
                    fprintf(['collecting inputs for ',input_stimtype,...
                        's ( SF=',num2str(pSF),' TF=',num2str(pTF),' ) neuron ',num2str(n),' collected ...\n']);
                    
                    toc
                    
                    % predict tuning curve and psth (testing) ------------------------------------------------------------------
                    
                    tic
                    
                    % compute predictions
                    sr=videoframerate;
                    countingwindowlim=[0,1];
                    alpha=alphaopt(i);
                    beta=betaopt(i);
                    [obs_fr,pred_fr,obs_count,pred_count,fr_time] =...
                        predict_neuron_responses(input_PSTHs,input_S,input_dirs,...
                        alpha,beta,input_F,sr,countingwindowlim);
                    % output progress message
                    fprintf(['filter response at ',input_stimtype,...
                        's ( SF=',num2str(pSF),' TF=',num2str(pTF),' ) neuron ',num2str(n),' computed ...\n']);
                    % get goodness of fit
                    [explvar_T(i,k),rmse_tuning,explvar_D(i,k),rmse_dynamics] = ...
                        get_neuron_responses_prediction_gof(obs_fr,obs_count,pred_fr,pred_count,DIR==pDIR);
                    % store results
                    obs_FR(:,:,i,k)=obs_fr;
                    obs_COUNT(:,i,k)=obs_count;
                    pred_FR(:,:,i,k)=pred_fr;
                    pred_COUNT(:,i,k)=pred_count;
                    FR_time=fr_time;
                    
                    % collect inputs for plotting
                    stimtoplot{k}=input_STIMULIs{pDIR==DIR};
                    rastertoplot{k}=input_RASTERs{pDIR==DIR};
                    
                    % clear heavy vars
                    clear input_F input_S input_STIMULIs
                    
                    toc
                    
                    % compute predicted pattern index (testing) ------------------------------------------------------------------
                    
                    %                 % get back observed index values
                    %                 PI_O(i)=PI(n,pSF==SF,pTF==TF);
                    %                 Zp_O(i)=Zp(n,pSF==SF,pTF==TF);
                    %                 Zc_O(i)=Zc(n,pSF==SF,pTF==TF);
                    %                 Rp_O(i)=Rp(n,pSF==SF,pTF==TF);
                    %                 Rc_O(i)=Rc(n,pSF==SF,pTF==TF);
                    %                 % reclassify
                    %                 if Zp_O(i)-max(Zc_O(i),0)>=1.28
                    %                     ctlabel_O(i)=2; % 2=pattern
                    %                 elseif Zc_O(i)-max(Zp_O(i),0)>=1.28
                    %                     ctlabel_O(i)=1; % 1=component
                    %                 else
                    %                     ctlabel_O(i)=0; % 0=unclassified
                    %                 end
                    
                    classlabels={'UNCL','COMP','PATT'};
                    if k==2
                        % get tuning curves
                        tuning_curve_grating_P=pred_COUNT(:,i,1);
                        tuning_curve_plaid_P=pred_COUNT(:,i,2);
                        % perform partial correlation analysis to get predicted index values
                        [ PI_P(i), ~, Zp_P(i), Zc_P(i), Rp_P(i), Rc_P(i), ~, ~, ~, ~ ] =...
                            get_pattern_index( tuning_curve_grating_P,tuning_curve_plaid_P );
                        % reclassify
                        if Zp_P(i)-max(Zc_P(i),0)>=1.28
                            ctlabel_P(i)=2; % 2=pattern
                        elseif Zc_P(i)-max(Zp_P(i),0)>=1.28
                            ctlabel_P(i)=1; % 1=component
                        else
                            ctlabel_P(i)=0; % 0=unclassified
                        end
                        % get tuning curves
                        tuning_curve_grating_O=obs_COUNT(:,i,1);
                        tuning_curve_plaid_O=obs_COUNT(:,i,2);
                        % perform partial correlation analysis to get predicted index values
                        [ PI_O(i), ~, Zp_O(i), Zc_O(i), Rp_O(i), Rc_O(i), ~, ~, ~, ~ ] =...
                            get_pattern_index( tuning_curve_grating_O,tuning_curve_plaid_O );
                        % reclassify
                        if Zp_O(i)-max(Zc_O(i),0)>=1.28
                            ctlabel_O(i)=2; % 2=pattern
                        elseif Zc_O(i)-max(Zp_O(i),0)>=1.28
                            ctlabel_O(i)=1; % 1=component
                        else
                            ctlabel_O(i)=0; % 0=unclassified
                        end
                    end
                    
                    % plot predictions ------------------------------------------------------------------
                    
                    tic
                    if k==2
                        for kk=1:2
                            
                            % get pref dir to plot
                            pDIR=pdirtoplot{kk};
                            pSF=psftoplot{kk};
                            pTF=ptftoplot{kk};
                            
                            % initialize figure
                            f1 = figure('units','normalized','outerposition',[0 0 1 1]);
                            % set psth subplot position ----------
                            sb1=subplot(666,666,1);
                            set(sb1,'Position',[.52,0.4,.5,.5]);
                            axis square
                            % get PSTHs to plot
                            psth_predicted=pred_FR(:,DIR==pDIR,i,kk);
                            psth_predicted=psth_predicted./(max(psth_predicted(:)));
                            psth_observed=obs_FR(:,DIR==pDIR,i,kk);
                            psth_observed=psth_observed./(max(psth_observed(:)));
                            % set color and tag to use
                            if ctlabel_O(i)==2
                                coltuse=[255,150,0]./255;
                            elseif ctlabel_O(i)==1
                                coltuse=[50,200,0]./255;
                            elseif ctlabel_O(i)==0
                                coltuse=[150,150,150]./255;
                            end
                            hold on;
                            % draw psth
                            plot(gca,FR_time,psth_observed,'-','Color',coltuse,'LineWidth',2.5);
                            plot_shaded_auc(gca,FR_time,psth_observed',0.15,coltuse)
                            plot(gca,FR_time,psth_predicted,':','Color',coltuse*0.5,'LineWidth',2.5);
                            % draw raster
                            raster=rastertoplot{kk};
                            spp=raster;
                            spk=0;
                            for kkk=1:numel(raster)
                                if not(isempty(raster{kkk}))
                                    plot(raster{kkk},1.05+0.05*kkk,'.k', 'MarkerSize',15)
                                    spk=spk+numel(raster{kkk});
                                else
                                end
                            end
                            xlim([-0.2,1.2]);
                            ylim([-0,4]);
                            plot([0,0],[0,5],'--k', 'LineWidth',2)
                            plot([1,1],[0,5],'--k', 'LineWidth',2)
                            tt=text(0.05,3.5,['DIR = ',num2str(pDIR),' d'],'FontSize',12);
                            ttt=text(0.05,3.25,['spike count = ',num2str(spk)],'FontSize',12);
                            set(gca,'FontSize',12);
                            tttt = text(0.05,3.00,['ev psth = ',num2str(explvar_D(i,kk),'%.2f')],'FontSize',12);
                            hlabelx=get(gca,'Xlabel');
                            set(hlabelx,'String','time (s)','FontSize',12,'color','k')
                            hlabely=get(gca,'Ylabel');
                            set(hlabely,'String','normalized firing rate','FontSize',12,'color','k')
                            % select bitcodes for the tuning curve
                            selected_idx = get_indexes_PN( pSF, pTF, pDIR, stimulustypes{kk} );
                            sessionname=[listSessions{1,M(n,1)},'_b',num2str(M(n,2))];
                            sname=strrep(sessionname,'_',' ');
                            title(['best ',stimulustypes{kk},' - neu ',num2str(n),' - session ',sname])
                            % set psth subplot position ----------
                            ppol=polaraxes('Position',[.02,0.4,.5,.5]);
                            % get tuning curve to plot
                            temp_pred_tc=squeeze(pred_COUNT(:,i,kk))'./max(squeeze(pred_COUNT(:,i,kk)));
                            temp_obs_tc=squeeze(obs_COUNT(:,i,kk))'./max(squeeze(obs_COUNT(:,i,kk)));
                            pred_tc=[temp_pred_tc,temp_pred_tc(1)];
                            obs_tc=[temp_obs_tc,temp_obs_tc(1)];
                            title(ppol,[stimulustypes{kk},' tuning curve - TF = ',num2str(pTF),' Hz - SF = ',num2str(pSF,'%.2f'),' cpd - avgCI = ',num2str(nanmean(cis_per_frame_to_use(:,i)))],'fontsize',12)
                            set(ppol,'fontsize',12);
                            hold on;
                            % draw plar plots
                            p1=polarplot(ppol,[deg2rad(DIR),2*pi],pred_tc,'-');
                            p2=polarplot(ppol,[deg2rad(DIR),2*pi],obs_tc,'-');
                            set(p1,'color',coltuse*0.5)
                            set(p2,'color',coltuse)
                            set(p1, 'linewidth', 3.5);
                            set(p2, 'linewidth', 3.5);
                            set(ppol,'fontsize',12)
                            tx0=text(ppol,deg2rad(45),1.5,['ev tuning = ',num2str(explvar_T(i,kk),'%.2f')],'fontsize',12);
                            tx1=text(ppol,deg2rad(40),1.5,['betaopt (npTF) = ',num2str(betaopt(i),'%.2f')],'fontsize',12);
                            % select best frame
                            [~,imaxfp]=max(r2s_per_frame_to_use(:,i));
                            %                 if k==2
                            tx2=text(ppol,deg2rad(30),1.5,['Zp obs = ',num2str(Zp_O(i),'%.01f')],'fontsize',12);
                            tx3=text(ppol,deg2rad(25),1.5,['Zp pred = ',num2str(Zp_P(i),'%.01f')],'fontsize',12);
                            tx4=text(ppol,deg2rad(20),1.5,['Zc obs = ',num2str(Zc_O(i),'%.01f')],'fontsize',12);
                            tx5=text(ppol,deg2rad(15),1.5,['Zc pred = ',num2str(Zc_P(i),'%.01f')],'fontsize',12);
                            tx6=text(ppol,deg2rad(10),1.5,['obs = ',classlabels{ctlabel_O(i)+1}],'fontsize',12);
                            tx7=text(ppol,deg2rad(5),1.5,['pred = ',classlabels{ctlabel_P(i)+1}],'fontsize',12);
                            %                 else
                            %                 end
                            % loop over frames to draw filter
                            for jj=1:size(selectedSTA,3)
                                sb3=subplot(666,666,1);
                                fram=selectedSTA_raw(:,:,jj);  % NB: added for plotting
                                %                                 if r2s_per_frame_to_use(jj,i)<=r2_th_for_pred_anlaysis
                                %                                     fram=zeros(size(fram));
                                %                                 end
                                fram=imresize(fram,3);
                                l1=imagesc(fram); colormap('gray');
                                % caxis([quantile(selectedSTA(:),0.01), quantile(selectedSTA(:),0.99)]); set(gca,'dataAspectRatio',[1 1 1]); axis off
                                caxis([-6,6]); set(gca,'dataAspectRatio',[1 1 1]); axis off   % NB: added for plotting
                                set(sb3,'Position',[.02+0.095*(jj-1),0.1,.09,.09]);
                                text(sb3,(2/3)*size(fram,1),60,['CI=',num2str(cis_per_frame_to_use(jj,i))]);  % NB: added for plotting
                                text(sb3,(2/3)*size(fram,1),70,['R2=',num2str(r2s_per_frame_to_use(jj,i))]);  % NB: added for plotting
                            end
                            hold on
                            % loop over frames to draw stimulus
                            pstimul=stimtoplot{kk};
                            for jj=1:size(selectedSTA,3)
                                sb4=subplot(666,666,1);
                                fram=pstimul(:,:,25+jj);
                                fram=imresize(fram,0.5);
                                imagesc(fram); colormap(gray); set(gca,'dataAspectRatio',[1 1 1]); axis off
                                set(sb4,'Position',[.02+0.095*(jj-1),0.2,.09,.09]);
                            end
                            % save
                            fname=['prediction ','(n_',num2str(n),' goodn_',num2str(i),') ',are,' ',stimulustypes{kk}];
                            fname=strrep(fname,'.','');
                            % NB: added for plotting
                            fname=[fname,'_PLOTTING'];
                            % save jpg(gcf,fname, 'epsc')
                            set(gcf, 'PaperPositionMode', 'auto')
                            saveas(gcf,fname, 'jpg')
                            % save eps NB: added for plotting
                            print(gcf,'-depsc','-painters',[fname,'.eps'])
                            
                            toc
                            
                        end
                        close all
                    end
                    
                end
                
                % NB: added for plotting (silenced saving below)
                % save results
                save([outfold,filesep,are,'_',cty,'_prediction_explained_variance_2022'],'pred_COUNT','obs_COUNT','explvar_T','explvar_D','explvar_sD','pDIRdiff','pORIdiff',...
                    'PI_O','Zp_O','Zc_O','Rp_O','Rc_O','PI_P','Zp_P','Zc_P','Rp_P','Rc_P','obs_OSI','obs_DSI','ctlabel_O','ctlabel_P','alphaopt','betaopt')
                %            end  % NB: added for plotting
            end
        end
    end
    
end

%% recollect goodness of gabor fit goodness, contrast factors and lobe counts

% load receptive field fitting data
rffit=load([outfold,filesep,'rf_fitting_results.mat']);
% get number of neurons per area and cell type
n_neurons=zeros(length(cell_types),length(target_areas));
% loop over areas
for target_area_idx=1:length(target_areas)%1:length(target_areas)
    % set area
    are=target_areas{target_area_idx};
    % loop over cell types
    for cell_type_idx=1:length(cell_types)
        n_neurons(cell_type_idx,target_area_idx)=size(rffit.fitted_rfs_r2_distribs{cell_type_idx,target_area_idx},2);
    end
end
% reorganize and collect lobe cout datastructure
best_lobe_number_distribs=cell(length(cell_types),length(target_areas));
% loop over cell types
for cell_type_idx=1:length(cell_types)
    % get indexes to grab per area
    indexes_to_grab{1}=1:n_neurons(cell_type_idx,1);
    indexes_to_grab{2}=n_neurons(cell_type_idx,1)+(1:n_neurons(cell_type_idx,2));
    indexes_to_grab{3}=sum(n_neurons(cell_type_idx,1:2))+(1:n_neurons(cell_type_idx,3));
    % loop over areas
    for target_area_idx=1:length(target_areas)
        best_lobe_number_distribs{cell_type_idx,target_area_idx}=...
            rffit.collapsed_best_lobe_number_distribs{cell_type_idx}(indexes_to_grab{target_area_idx});
    end
end
% collect other datastructures
fitted_rfs_r2_distribs=rffit.fitted_rfs_r2_distribs;
contrast_index_distribs=rffit.contrast_index_distribs;
% save addendum data structute
for target_area_idx=1:length(target_areas)%1:length(target_areas)
    % set area
    are=target_areas{target_area_idx};
    % loop over cell types
    for cell_type_idx=1:length(cell_types)
        tic
        % set cell type
        cty=cell_types{cell_type_idx};
        % load neuron selection results
        sel_res_file=['D:\Backups\Personal_bk\PN_acute_analysis\processed_data',filesep,are,'_selected_neurons',filesep,'selected_population_',are,'.mat'];
        load(sel_res_file)
        % load fitted sta analysis results
        sta_res_file=[outfold,filesep,'rf_fitting_results.mat'];
        load(sta_res_file)
        % select current category neurons with good enough Gabor fitted rf
        current_r2s_per_frame=fitted_rfs_r2_distribs{cell_type_idx,target_area_idx};
        valid_neurons_idx=find(max(current_r2s_per_frame,[],1)>r2_th_for_pred_anlaysis);
        current_selected_neu_idx=selected_cells_idx_distribs{cell_type_idx,target_area_idx}(valid_neurons_idx);
        selectedsi_to_use=selectedsi(current_selected_neu_idx); % inex within the ogiginal good .. structure
        % recover original STA for display
        sta_res_file=[data_folder,filesep,'RFs_new_datastructure_',are];
        load(sta_res_file)
        cell_types_codes=[1,2,0];
        current_cells_idx=find(goomlabel==cell_types_codes(cell_type_idx));
        cis_per_frame_to_use=contrast(:,current_cells_idx(valid_neurons_idx));
        % get old best condition Zp and Zc
        Zc_O_oldbest=goozc(current_selected_neu_idx);
        Zp_O_oldbest=goozp(current_selected_neu_idx);
        % recover old pattern and component index for comparison
        PI_O_old=zeros(numel(selectedsi_to_use),1);
        Zp_O_old=zeros(numel(selectedsi_to_use),1);
        Zc_O_old=zeros(numel(selectedsi_to_use),1);
        Rp_O_old=zeros(numel(selectedsi_to_use),1);
        Rc_O_old=zeros(numel(selectedsi_to_use),1);
        ns_old=zeros(numel(selectedsi_to_use),1);
        pTFs_old=zeros(numel(selectedsi_to_use),1);
        pSFs_old=zeros(numel(selectedsi_to_use),1);
        for iii=1:numel(selectedsi_to_use)
            pSF=SF(1); %SF(1);
            pTF=TF(selectedsi_to_use{iii}(3));
            n=selectedsi_to_use{iii}(1);
            PI_O_old(iii)=PI(n,pSF==SF,pTF==TF);
            Zp_O_old(iii)=Zp(n,pSF==SF,pTF==TF);
            Zc_O_old(iii)=Zc(n,pSF==SF,pTF==TF);
            Rp_O_old(iii)=Rp(n,pSF==SF,pTF==TF);
            Rc_O_old(iii)=Rc(n,pSF==SF,pTF==TF);
            ns_old(iii)=n;
            pTFs_old(iii)=TF(selectedsi_to_use{iii}(3));
            pSFs_old(iii)=SF(selectedsi_to_use{iii}(2));
        end
        % get current gabor fit goodness, contrast factors, relative areas and lobe counts
        lobecount=best_lobe_number_distribs{cell_type_idx,target_area_idx};
        gaborfitr2=fitted_rfs_r2_distribs{cell_type_idx,target_area_idx}';
        contrastindex=cis_per_frame_to_use';%contrast_index_distribs{cell_type_idx,target_area_idx};
        % save addendum file to prediction one containing pld indexes and rf quantification
        save([outfold,filesep,are,'_',cty,'_prediction_explained_variance_addendum_2022'],...
            'PI_O_old',...
            'Zp_O_old',...
            'Zc_O_old',...
            'Zp_O_oldbest',...
            'Zc_O_oldbest',...
            'Rp_O_old',...
            'Rc_O_old',...
            'ns_old',...
            'pTFs_old',...
            'pSFs_old',...
            'lobecount',...
            'gaborfitr2',...
            'contrastindex'...
            );
        % TODO: add storage of absolute neuron number n
        % selectedsi_to_use(1) to see wich ones are selected at the end and
        % check which single neuron examples are usable
        toc
    end
end

%% collect tuning prediction results

%NB: this is also where cells gets reclassified according to the new Zp Zc

% outfold=[outfold,'\28_09_2022_bk'];
% save('test.mat',...
%     'Zc_O_distribs','Zc_O_distribs');

% initialize storage
explvart_distribs=cell(length(cell_types),3);
explvard_distribs=cell(length(cell_types),3);
Zc_O_distribs=cell(length(cell_types),3);
Zp_O_distribs=cell(length(cell_types),3);
Zc_O_old_distribs=cell(length(cell_types),3);
Zp_O_old_distribs=cell(length(cell_types),3);
Zc_O_oldbest_distribs=cell(length(cell_types),3);
Zp_O_oldbest_distribs=cell(length(cell_types),3);
Zc_P_distribs=cell(length(cell_types),3);
Zp_P_distribs=cell(length(cell_types),3);
PI_O_distribs=cell(length(cell_types),3);
PI_O_old_distribs=cell(length(cell_types),3);
PI_P_distribs=cell(length(cell_types),3);
lobecount_distribs=cell(length(cell_types),3);
ns_distribs=cell(length(cell_types),3);
pSF_distribs=cell(length(cell_types),3);
pTF_distribs=cell(length(cell_types),3);
gaborfitr2_distribs=cell(length(cell_types),3);
contrastindex_distribs=cell(length(cell_types),3);
Explvar_Tg_distribs=cell(length(cell_types),3);
Explvar_Dg_distribs=cell(length(cell_types),3);
Explvar_Tp_distribs=cell(length(cell_types),3);
Explvar_Dp_distribs=cell(length(cell_types),3);
pred_COUNT_distribs=cell(length(cell_types),3);
obs_COUNT_distribs=cell(length(cell_types),3);
obs_OSI_distribs=cell(length(cell_types),3);
obs_DSI_distribs=cell(length(cell_types),3);
ctlabel_O_distribs=cell(length(cell_types),3);
ctlabel_P_distribs=cell(length(cell_types),3);
alphaopt_distribs=cell(length(cell_types),3);
betaopt_distribs=cell(length(cell_types),3);
% loop over areas
target_areas={'V1','LM','RL'};
cell_types={'component','pattern','unclassified'};
cell_types_codes=[1,2,0];
for target_area_idx=1:length(target_areas)
    % set area
    are=target_areas{target_area_idx};
    % loop over original cell types
    for cell_type_idx=1:length(cell_types)
        % set original cell type
        cty=cell_types{cell_type_idx};
        % load results
        rfpred=load([outfold,filesep,are,'_',cty,'_prediction_explained_variance_2022.mat']);
        rffit=load([outfold,filesep,are,'_',cty,'_prediction_explained_variance_addendum_2022.mat']);
        % decide wich classification to use to resort cells
        switch labelstouse
            case 'oldbest'
                bool_above_zp=(rffit.Zp_O_oldbest-max(rffit.Zc_O_oldbest,zeros(size(rffit.Zc_O_oldbest))))>=1.28;
                bool_above_zc=(rffit.Zc_O_oldbest-max(rffit.Zp_O_oldbest,zeros(size(rffit.Zp_O_oldbest))))>=1.28;
                ctindex=zeros(size(bool_above_zp));
                ctindex(bool_above_zc)=1;
                ctindex(bool_above_zp)=2;
            case 'old'
                bool_above_zp=(rffit.Zp_O_old-max(rffit.Zc_O_old,zeros(size(rffit.Zc_O_old))))>=1.28;
                bool_above_zc=(rffit.Zc_O_old-max(rffit.Zp_O_old,zeros(size(rffit.Zp_O_old))))>=1.28;
                ctindex=zeros(size(bool_above_zp));
                ctindex(bool_above_zc)=1;
                ctindex(bool_above_zp)=2;
            case 'new'
                ctindex=rfpred.ctlabel_O;
            case 'oldbest_intersection'
                bool_above_zp=(rffit.Zp_O_oldbest-max(rffit.Zc_O_oldbest,zeros(size(rffit.Zc_O_oldbest))))>=1.28;
                bool_above_zc=(rffit.Zc_O_oldbest-max(rffit.Zp_O_oldbest,zeros(size(rffit.Zp_O_oldbest))))>=1.28;
                ctindex=zeros(size(bool_above_zp));
                ctindex(bool_above_zc)=1;
                ctindex(bool_above_zp)=2;
                ctindex(rfpred.ctlabel_O~=ctindex)=4;
        end
        % loop over cell classes
        for cell_type_idx_bis=1:length(cell_types)
            curr_cell_idx=find(ctindex==cell_types_codes(cell_type_idx_bis));
            % stack in variables
            explvart_distribs{cell_type_idx_bis,target_area_idx}=cat(1,explvart_distribs{cell_type_idx_bis,target_area_idx},rfpred.explvar_T(curr_cell_idx,:));
            explvard_distribs{cell_type_idx_bis,target_area_idx}=cat(1,explvard_distribs{cell_type_idx_bis,target_area_idx},rfpred.explvar_D(curr_cell_idx,:));
            Zc_O_distribs{cell_type_idx_bis,target_area_idx}=cat(1,Zc_O_distribs{cell_type_idx_bis,target_area_idx},rfpred.Zc_O(curr_cell_idx));
            Zp_O_distribs{cell_type_idx_bis,target_area_idx}=cat(1,Zp_O_distribs{cell_type_idx_bis,target_area_idx},rfpred.Zp_O(curr_cell_idx));
            Zc_O_old_distribs{cell_type_idx_bis,target_area_idx}=cat(1,Zc_O_old_distribs{cell_type_idx_bis,target_area_idx},rffit.Zc_O_old(curr_cell_idx));
            Zp_O_old_distribs{cell_type_idx_bis,target_area_idx}=cat(1,Zp_O_old_distribs{cell_type_idx_bis,target_area_idx},rffit.Zp_O_old(curr_cell_idx));
            Zc_O_oldbest_distribs{cell_type_idx_bis,target_area_idx}=cat(1,Zc_O_oldbest_distribs{cell_type_idx_bis,target_area_idx},rffit.Zc_O_oldbest(curr_cell_idx));
            Zp_O_oldbest_distribs{cell_type_idx_bis,target_area_idx}=cat(1,Zp_O_oldbest_distribs{cell_type_idx_bis,target_area_idx},rffit.Zp_O_oldbest(curr_cell_idx));
            Zc_P_distribs{cell_type_idx_bis,target_area_idx}=cat(1,Zc_P_distribs{cell_type_idx_bis,target_area_idx},rfpred.Zc_P(curr_cell_idx));
            Zp_P_distribs{cell_type_idx_bis,target_area_idx}=cat(1,Zp_P_distribs{cell_type_idx_bis,target_area_idx},rfpred.Zp_P(curr_cell_idx));
            PI_O_distribs{cell_type_idx_bis,target_area_idx}=cat(1,PI_O_distribs{cell_type_idx_bis,target_area_idx},rfpred.PI_O(curr_cell_idx));
            PI_O_old_distribs{cell_type_idx_bis,target_area_idx}=cat(1,PI_O_old_distribs{cell_type_idx_bis,target_area_idx},rffit.PI_O_old(curr_cell_idx));
            PI_P_distribs{cell_type_idx_bis,target_area_idx}=cat(1,PI_P_distribs{cell_type_idx_bis,target_area_idx},rfpred.PI_P(curr_cell_idx));
            Explvar_Tg_distribs{cell_type_idx_bis,target_area_idx}=cat(1,Explvar_Tg_distribs{cell_type_idx_bis,target_area_idx},rfpred.explvar_T(curr_cell_idx,1));
            Explvar_Dg_distribs{cell_type_idx_bis,target_area_idx}=cat(1,Explvar_Dg_distribs{cell_type_idx_bis,target_area_idx},rfpred.explvar_D(curr_cell_idx,1));
            Explvar_Tp_distribs{cell_type_idx_bis,target_area_idx}=cat(1,Explvar_Tp_distribs{cell_type_idx_bis,target_area_idx},rfpred.explvar_T(curr_cell_idx,2));
            Explvar_Dp_distribs{cell_type_idx_bis,target_area_idx}=cat(1,Explvar_Dp_distribs{cell_type_idx_bis,target_area_idx},rfpred.explvar_D(curr_cell_idx,2));
            lobecount_distribs{cell_type_idx_bis,target_area_idx}=cat(1,lobecount_distribs{cell_type_idx_bis,target_area_idx},rffit.lobecount(curr_cell_idx));
            ns_distribs{cell_type_idx_bis,target_area_idx}=cat(1,ns_distribs{cell_type_idx_bis,target_area_idx},rffit.ns_old(curr_cell_idx));
            pSF_distribs{cell_type_idx_bis,target_area_idx}=cat(1,pSF_distribs{cell_type_idx_bis,target_area_idx},rffit.pSFs_old(curr_cell_idx));
            pTF_distribs{cell_type_idx_bis,target_area_idx}=cat(1,pTF_distribs{cell_type_idx_bis,target_area_idx},rffit.pTFs_old(curr_cell_idx));
            gaborfitr2_distribs{cell_type_idx_bis,target_area_idx}=cat(1,gaborfitr2_distribs{cell_type_idx_bis,target_area_idx},rffit.gaborfitr2(curr_cell_idx,:));
            contrastindex_distribs{cell_type_idx_bis,target_area_idx}=cat(1,contrastindex_distribs{cell_type_idx_bis,target_area_idx},rffit.contrastindex(curr_cell_idx,:));
            temp=rfpred.pred_COUNT(:,curr_cell_idx,:);
            if not(isempty(temp))
                pred_COUNT_distribs{cell_type_idx_bis,target_area_idx}=cat(2,pred_COUNT_distribs{cell_type_idx_bis,target_area_idx},temp);
            end
            temp=rfpred.obs_COUNT(:,curr_cell_idx,:);
            if not(isempty(temp))
                obs_COUNT_distribs{cell_type_idx_bis,target_area_idx}=cat(2,obs_COUNT_distribs{cell_type_idx_bis,target_area_idx},temp);
            end
            obs_OSI_distribs{cell_type_idx_bis,target_area_idx}=cat(1,obs_OSI_distribs{cell_type_idx_bis,target_area_idx},rfpred.obs_OSI(curr_cell_idx));
            obs_DSI_distribs{cell_type_idx_bis,target_area_idx}=cat(1,obs_DSI_distribs{cell_type_idx_bis,target_area_idx},rfpred.obs_DSI(curr_cell_idx));
        end
    end
end
% initialize storage
collapsed_explvartG_distribs=cell(length(cell_types),1);
collapsed_explvardG_distribs=cell(length(cell_types),1);
collapsed_explvartP_distribs=cell(length(cell_types),1);
collapsed_explvardP_distribs=cell(length(cell_types),1);
collapsed_explvart_distribs=cell(length(cell_types),1);
collapsed_explvard_distribs=cell(length(cell_types),1);
collapsed_Zc_O_distribs=cell(length(cell_types),1);
collapsed_Zp_O_distribs=cell(length(cell_types),1);
collapsed_Zc_O_old_distribs=cell(length(cell_types),1);
collapsed_Zp_O_old_distribs=cell(length(cell_types),1);
collapsed_Zc_O_oldbest_distribs=cell(length(cell_types),1);
collapsed_Zp_O_oldbest_distribs=cell(length(cell_types),1);
collapsed_Zc_P_distribs=cell(length(cell_types),1);
collapsed_Zp_P_distribs=cell(length(cell_types),1);
collapsed_PI_O_distribs=cell(length(cell_types),1);
collapsed_PI_O_old_distribs=cell(length(cell_types),1);
collapsed_PI_P_distribs=cell(length(cell_types),1);
collapsed_Explvar_Tg_distribs=cell(length(cell_types),1);
collapsed_Explvar_Dg_distribs=cell(length(cell_types),1);
collapsed_Explvar_Tp_distribs=cell(length(cell_types),1);
collapsed_Explvar_Dp_distribs=cell(length(cell_types),1);
collapsed_pred_COUNT_distribs=cell(length(cell_types),1);
collapsed_obs_COUNT_distribs=cell(length(cell_types),1);
collapsed_obs_OSI_distribs=cell(length(cell_types),1);
collapsed_obs_DSI_distribs=cell(length(cell_types),1);
collapsed_lobecount_distribs=cell(length(cell_types),1);
collapsed_ns_distribs=cell(length(cell_types),1);
collapsed_pSF_distribs=cell(length(cell_types),1);
collapsed_pTF_distribs=cell(length(cell_types),1);
collapsed_gaborfitr2_distribs=cell(length(cell_types),1);
collapsed_contrastindex_distribs=cell(length(cell_types),1);
% loop over cell types
for cell_type_idx=1:length(cell_types)
    % loop over areas
    for target_area_idx=1:length(target_areas)
        collapsed_explvartG_distribs{cell_type_idx}=cat(1,collapsed_explvartG_distribs{cell_type_idx},explvart_distribs{cell_type_idx,target_area_idx}(:,1));
        collapsed_explvardG_distribs{cell_type_idx}=cat(1,collapsed_explvardG_distribs{cell_type_idx},explvard_distribs{cell_type_idx,target_area_idx}(:,1));
        collapsed_explvartP_distribs{cell_type_idx}=cat(1,collapsed_explvartP_distribs{cell_type_idx},explvart_distribs{cell_type_idx,target_area_idx}(:,2));
        collapsed_explvardP_distribs{cell_type_idx}=cat(1,collapsed_explvardP_distribs{cell_type_idx},explvard_distribs{cell_type_idx,target_area_idx}(:,2));
        collapsed_explvart_distribs{cell_type_idx}=cat(1,collapsed_explvart_distribs{cell_type_idx},explvart_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_explvard_distribs{cell_type_idx}=cat(1,collapsed_explvard_distribs{cell_type_idx},explvard_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_Zc_O_distribs{cell_type_idx}=cat(1,collapsed_Zc_O_distribs{cell_type_idx},Zc_O_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_Zp_O_distribs{cell_type_idx}=cat(1,collapsed_Zp_O_distribs{cell_type_idx},Zp_O_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_Zc_O_old_distribs{cell_type_idx}=cat(1,collapsed_Zc_O_old_distribs{cell_type_idx},Zc_O_old_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_Zp_O_old_distribs{cell_type_idx}=cat(1,collapsed_Zp_O_old_distribs{cell_type_idx},Zp_O_old_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_Zc_O_oldbest_distribs{cell_type_idx}=cat(1,collapsed_Zc_O_oldbest_distribs{cell_type_idx},Zc_O_oldbest_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_Zp_O_oldbest_distribs{cell_type_idx}=cat(1,collapsed_Zp_O_oldbest_distribs{cell_type_idx},Zp_O_oldbest_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_Zc_P_distribs{cell_type_idx}=cat(1,collapsed_Zc_P_distribs{cell_type_idx},Zc_P_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_Zp_P_distribs{cell_type_idx}=cat(1,collapsed_Zp_P_distribs{cell_type_idx},Zp_P_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_PI_O_distribs{cell_type_idx}=cat(1,collapsed_PI_O_distribs{cell_type_idx},PI_O_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_PI_O_old_distribs{cell_type_idx}=cat(1,collapsed_PI_O_old_distribs{cell_type_idx},PI_O_old_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_PI_P_distribs{cell_type_idx}=cat(1,collapsed_PI_P_distribs{cell_type_idx},PI_P_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_Explvar_Tg_distribs{cell_type_idx}=cat(1,collapsed_Explvar_Tg_distribs{cell_type_idx},Explvar_Tg_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_Explvar_Dg_distribs{cell_type_idx}=cat(1,collapsed_Explvar_Dg_distribs{cell_type_idx},Explvar_Dg_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_Explvar_Tp_distribs{cell_type_idx}=cat(1,collapsed_Explvar_Tp_distribs{cell_type_idx},Explvar_Tp_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_Explvar_Dp_distribs{cell_type_idx}=cat(1,collapsed_Explvar_Dp_distribs{cell_type_idx},Explvar_Dp_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_pred_COUNT_distribs{cell_type_idx}=cat(2,collapsed_pred_COUNT_distribs{cell_type_idx},pred_COUNT_distribs{cell_type_idx,target_area_idx});
        collapsed_obs_COUNT_distribs{cell_type_idx}=cat(2,collapsed_obs_COUNT_distribs{cell_type_idx},obs_COUNT_distribs{cell_type_idx,target_area_idx});
        collapsed_obs_OSI_distribs{cell_type_idx}=cat(1,collapsed_obs_OSI_distribs{cell_type_idx},obs_OSI_distribs{cell_type_idx,target_area_idx});
        collapsed_obs_DSI_distribs{cell_type_idx}=cat(1,collapsed_obs_DSI_distribs{cell_type_idx},obs_DSI_distribs{cell_type_idx,target_area_idx});
        collapsed_lobecount_distribs{cell_type_idx}=cat(1,collapsed_lobecount_distribs{cell_type_idx},lobecount_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_ns_distribs{cell_type_idx}=cat(1,collapsed_ns_distribs{cell_type_idx},ns_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_pSF_distribs{cell_type_idx}=cat(1,collapsed_pSF_distribs{cell_type_idx},pSF_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_pTF_distribs{cell_type_idx}=cat(1,collapsed_pTF_distribs{cell_type_idx},pTF_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_gaborfitr2_distribs{cell_type_idx}=cat(1,collapsed_gaborfitr2_distribs{cell_type_idx},gaborfitr2_distribs{cell_type_idx,target_area_idx}(:));
        collapsed_contrastindex_distribs{cell_type_idx}=cat(1,collapsed_contrastindex_distribs{cell_type_idx},contrastindex_distribs{cell_type_idx,target_area_idx}(:));
    end
end

%% plot reclassification diagnostics

% compute overlap between classifications
Zp_oldbest=[collapsed_Zp_O_oldbest_distribs{1};collapsed_Zp_O_oldbest_distribs{2};collapsed_Zp_O_oldbest_distribs{3}];
Zc_oldbest=[collapsed_Zc_O_oldbest_distribs{1};collapsed_Zc_O_oldbest_distribs{2};collapsed_Zc_O_oldbest_distribs{3}];
Zp_new=[collapsed_Zp_O_distribs{1};collapsed_Zp_O_distribs{2};collapsed_Zp_O_distribs{3}];
Zc_new=[collapsed_Zc_O_distribs{1};collapsed_Zc_O_distribs{2};collapsed_Zc_O_distribs{3}];
% get all class label with old best condition Zc and Zp
bool_above_zp=(Zp_oldbest-max(Zc_oldbest,zeros(size(Zc_oldbest))))>=1.28;
bool_above_zc=(Zc_oldbest-max(Zp_oldbest,zeros(size(Zp_oldbest))))>=1.28;
ctindexes_oldbest=zeros(size(bool_above_zp));
ctindexes_oldbest(bool_above_zc)=1;
ctindexes_oldbest(bool_above_zp)=2;
% get all class label with new Zc and Zp
bool_above_zp=(Zp_new-max(Zc_new,zeros(size(Zc_new))))>=1.28;
bool_above_zc=(Zc_new-max(Zp_new,zeros(size(Zp_new))))>=1.28;
ctindexes_new=zeros(size(bool_above_zp));
ctindexes_new(bool_above_zc)=1;
ctindexes_new(bool_above_zp)=2;
% compute intersection and overlaps
compidx_oldbest=find(ctindexes_oldbest==1);
compidx_new=find(ctindexes_new==1);
compidx_intersection=intersect(compidx_oldbest,compidx_new);
pattidx_oldbest=find(ctindexes_oldbest==2);
pattidx_new=find(ctindexes_new==2);
pattidx_both=intersect(pattidx_oldbest,pattidx_new);
pattidx_oldonly=setdiff(pattidx_oldbest,pattidx_new);
pattidx_newonly=setdiff(pattidx_new,pattidx_oldbest);
compidx_both=intersect(compidx_oldbest,compidx_new);
compidx_oldonly=setdiff(compidx_oldbest,compidx_new);
compidx_newonly=setdiff(compidx_new,compidx_oldbest);
% initialize figure
compcol=[50,200,0]./255;
pattcol=[255,150,0]./255;
fighanddiag=figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,2)
hold on;
oldcompidx=(ctindexes_oldbest==1);
oldpattidx=(ctindexes_oldbest==2);
newcompidx=(ctindexes_new==1);
newpattidx=(ctindexes_new==2);
scatter(Zp_oldbest(oldcompidx),Zp_new(oldcompidx),100,...
    'markerfacecolor',compcol,'markeredgecolor',compcol)
scatter(Zp_oldbest(oldpattidx),Zp_new(oldpattidx),100,...
    'markerfacecolor',pattcol,'markeredgecolor',pattcol)
scatter(Zp_oldbest(newcompidx),Zp_new(newcompidx),50,...
    'markerfacecolor',compcol/2,'markeredgecolor',compcol/2)
scatter(Zp_oldbest(newpattidx),Zp_new(newpattidx),50,...
    'markerfacecolor',pattcol/2,'markeredgecolor',pattcol/2)
[corrval_Zp,corrp_Zp] = corr(Zp_oldbest,Zp_new);
ylimtouse=[-4,6];
xlimtouse=[-4,6];
ylim(gca,ylimtouse)
xlim(gca,ylimtouse)
plot(gca,xlimtouse,ylimtouse,'--','linewidth',2,'color',[0.5,0.5,0.5])
plot(gca,xlimtouse-1.26,ylimtouse,':','linewidth',2,'color',[0.5,0.5,0.5])
plot(gca,xlimtouse+1.26,ylimtouse,':','linewidth',2,'color',[0.5,0.5,0.5])
text(xlimtouse(1)+0.5,0.85*xlimtouse(2),['old vs new corr = ',num2str(corrval_Zp)],'fontsize',12)
title(gca,'Zp - old vs. new')
xlabel(gca,'old Zp (best SF)')
ylabel(gca,'new Zp (low SF)')
set(gca,'fontsize',12)
axis square
subplot(1,3,1)
hold on;
oldcompidx=(ctindexes_oldbest==1);
oldpattidx=(ctindexes_oldbest==2);
newcompidx=(ctindexes_new==1);
newpattidx=(ctindexes_new==2);
scatter(Zc_oldbest(oldcompidx),Zc_new(oldcompidx),100,...
    'markerfacecolor',compcol,'markeredgecolor',compcol)
scatter(Zc_oldbest(oldpattidx),Zc_new(oldpattidx),100,...
    'markerfacecolor',pattcol,'markeredgecolor',pattcol)
scatter(Zc_oldbest(newcompidx),Zc_new(newcompidx),50,...
    'markerfacecolor',compcol/2,'markeredgecolor',compcol/2)
scatter(Zc_oldbest(newpattidx),Zc_new(newpattidx),50,...
    'markerfacecolor',pattcol/2,'markeredgecolor',pattcol/2)
[corrval_Zc,corrp_Zc] = corr(Zc_oldbest,Zc_new);
ylimtouse=[-4,6];
xlimtouse=[-4,6];
ylim(gca,ylimtouse)
xlim(gca,ylimtouse)
plot(gca,xlimtouse,ylimtouse,'--','linewidth',2,'color',[0.5,0.5,0.5])
plot(gca,xlimtouse-1.26,ylimtouse,':','linewidth',2,'color',[0.5,0.5,0.5])
plot(gca,xlimtouse+1.26,ylimtouse,':','linewidth',2,'color',[0.5,0.5,0.5])
text(xlimtouse(1)+0.5,0.85*xlimtouse(2),['old vs new corr = ',num2str(corrval_Zc)],'fontsize',12)
title(gca,'Zc - old vs. new')
xlabel(gca,'old Zc (best SF)')
ylabel(gca,'new Zc (low SF)')
set(gca,'fontsize',12)
axis square
subplot(1,3,3)
hold on;
pattsets=[numel(pattidx_both);numel(pattidx_oldonly);numel(pattidx_newonly)];
compsets=[numel(compidx_both);numel(compidx_oldonly);numel(compidx_newonly)];
bar([1,4,6],compsets./sum(compsets),'facecolor',compcol,'edgecolor',compcol,'barwidth',0.3)
bar([2,5,7],pattsets./sum(pattsets),'facecolor',pattcol,'edgecolor',pattcol,'barwidth',0.3)
xtickstouse=[[1,4,6],[2,5,7]];
xticks(sort(xtickstouse));
xzickslabelsss={'comp common','patt common','comp oldonly','patt oldonly','comp newonly','patt newonly'};
xticklabels(xzickslabelsss);
title(gca,['fraction of overlap/difference - ( patt c = ',...
    num2str(numel(pattidx_both)),' comp c = ',num2str(numel(compidx_both)),' )'])
xtickangle(45)
ylabel(gca,'fraction')
set(gca,'fontsize',12)
axis square
saveas(fighanddiag,[outfold,filesep,'reclassification_diagnostics_old_new'],'jpg')
print(fighanddiag,'-depsc','-painters',[[outfold,filesep,'reclassification_diagnostics_old_new'],'.eps'])

% set single neuron examples (manually choosen)
% examples_comp_n=[44,403,625];
% examples_patt_n=[137,396,427,438,532,547];
examples_comp_n=[13,15,48,162,175,87,376,403,436,587,603,625,809,1097,1180];
examples_patt_n=[117,419,530];
% verify that choosen single neuron examples are in the pool
examples_comp_n_present=intersect(examples_comp_n,collapsed_ns_distribs{1});
examples_patt_n_present=intersect(examples_patt_n,collapsed_ns_distribs{2});
examples_idx_present{1}=zeros(1,numel(examples_comp_n_present));
for ii=1:numel(examples_comp_n_present)
    examples_idx_present{1}(ii)=find(examples_comp_n_present(ii)==collapsed_ns_distribs{1});
end
examples_idx_present{2}=zeros(1,numel(examples_patt_n_present));
for ii=1:numel(examples_patt_n_present)
    examples_idx_present{2}(ii)=find(examples_patt_n_present(ii)==collapsed_ns_distribs{2});
end

%% plot shape features distributions violins

clear ylabellist yimtouselist ks_ban
distrtoplotlist_orig{1}=collapsed_contrastindex_distribs;
distrtoplotlist_orig{2}=collapsed_gaborfitr2_distribs;
distrtoplotlist_orig{3}=collapsed_lobecount_distribs;
distrtoplotlist{1}=cell(size(collapsed_contrastindex_distribs));
distrtoplotlist{2}=cell(size(collapsed_contrastindex_distribs));
distrtoplotlist{3}=cell(size(collapsed_contrastindex_distribs));
% filter according to selectivity
th_val=0;
for ii=1:numel(distrtoplotlist_orig{1})
    temp=repmat(collapsed_obs_DSI_distribs{ii},[1,10])';
    distribtouse_filtvals=temp(:);
    distrtoplotlist{1}{ii}=distrtoplotlist_orig{1}{ii}(distribtouse_filtvals>=th_val);
end
for ii=1:numel(distrtoplotlist_orig{2})
    temp=repmat(collapsed_obs_DSI_distribs{ii},[1,10])';
    distribtouse_filtvals=temp(:);
    distrtoplotlist{2}{ii}=distrtoplotlist_orig{2}{ii}(distribtouse_filtvals>=th_val);
end
for ii=1:numel(distrtoplotlist_orig{3})
    temp=repmat(collapsed_obs_DSI_distribs{ii},[1,1])';
    distribtouse_filtvals=temp(:);
    distrtoplotlist{3}{ii}=distrtoplotlist_orig{3}{ii}(distribtouse_filtvals>=th_val);
end
ylabellist{1}='contrast factor';
ylabellist{2}='Gabor r2';
yimtouselist{1}=[2.5,20];
yimtouselist{2}=[0.025,0.75];
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
    inputaxh = plot_violinplot_PN_new(inputadata,inputpars); %#ok<NASGU>
    pvalw = ranksum(distribtouse{1},distribtouse{2});
    [junk,pvalt] = ttest2(distribtouse{1},distribtouse{2});
    text(-0.2,0.95*max(yimtouselist{jj}),['median diff p = ',num2str(pvalw)],'fontsize',12);
    text(-0.2,0.92*max(yimtouselist{jj}),['mean diff p = ',num2str(pvalt)],'fontsize',12);
    set(gca,'fontsize',12)
    axis square
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
[~,p_chitest] = chiSquareTest([N_comp_lobe_n;N_patt_lobe_n]);
text(20,0.92,['chi square p = ',num2str(p_chitest)],'fontsize',12);
xticks(10*X+1.5)
xticklabels(X)
xlabel('# of lobes')
ylabel('fraction of cells')
set(gca,'fontsize',12)
axis square
sgtitle(['RF shape analysis - all areas'])
saveas(fighand2,[outfold,filesep,'RF_shape_analysis_all_areas'],'jpg')
print(fighand2,'-depsc','-painters',[[outfold,filesep,'RF_shape_analysis_all_areas'],'.eps'])

%% visualize tuning prediction comparison grating, plaid, pattern, component - r2

distribtouse{1}=collapsed_Explvar_Tg_distribs{1};
distribtouse{2}=collapsed_Explvar_Tp_distribs{1};
distribtouse{3}=collapsed_Explvar_Tg_distribs{2};
distribtouse{4}=collapsed_Explvar_Tp_distribs{2};
[~,p_ttestpr2_comp] = ttest(collapsed_Explvar_Tg_distribs{1},collapsed_Explvar_Tp_distribs{1});
[~,p_ttestpr2_patt] = ttest(collapsed_Explvar_Tg_distribs{2},collapsed_Explvar_Tp_distribs{2});
[~,p_ttestr2_cpg] = ttest2(collapsed_Explvar_Tg_distribs{1},collapsed_Explvar_Tg_distribs{2});
[~,p_ttestgr2_cpp] = ttest2(collapsed_Explvar_Tp_distribs{1},collapsed_Explvar_Tp_distribs{2});
Explvar_T_comp_diff=-(collapsed_Explvar_Tg_distribs{1}-collapsed_Explvar_Tp_distribs{1});
Explvar_T_patt_diff=-(collapsed_Explvar_Tg_distribs{2}-collapsed_Explvar_Tp_distribs{2});
ylabellist='tuning explained variance';
% yimtouselist=[-0.1,0.8];
% ks_ban=0.05;
yimtouselist=[-0.1,1.1];
ks_ban=0.1;%0.10;
[hh,pp]=ttest2(Explvar_T_comp_diff,Explvar_T_patt_diff);
% % % % perform anova
% % % obsanova=[distribtouse{1}',...
% % %     distribtouse{2}',...
% % %     distribtouse{3}',....
% % %     distribtouse{4}'];
% % % groupsanova=[repmat({'1'},[1,numel(distribtouse{1})]),...
% % %     repmat({'1'},[1,numel(distribtouse{2})]),...
% % %     repmat({'2'},[1,numel(distribtouse{3})]),...
% % %     repmat({'2'},[1,numel(distribtouse{4})])];
% % % panovar2tuning = anova1(obsanova,groupsanova);
% % % [hhhh,pppp]=ttest2([distribtouse{1};distribtouse{2}],[distribtouse{2};distribtouse{3}]);

% initialize figure
fighand1=figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
inputpars.inputaxh=gca;
hold(inputpars.inputaxh,'on')
% set settings for violin distribution plotting
inputpars.boxplotwidth=0.4;
inputpars.boxplotlinewidth=2;
inputpars.densityplotwidth=0.4;
inputpars.yimtouse=yimtouselist;
% inputpars.yimtouse=[0,8];
inputpars.scatterjitter=inputpars.boxplotlinewidth*0.1;
inputpars.scatteralpha=0.15;
inputpars.scattersize=50;
inputpars.distralpha=0.5;
inputpars.xlabelstring=[];
inputpars.ylabelstring=ylabellist;
inputpars.titlestring=[ylabellist,' p vs. g ( comp = ',...
    num2str(numel(distribtouse{1}),'%.0f'),...
    ' - patt  = ',num2str(numel(distribtouse{3}),'%.0f'),' )'];
inputpars.boolscatteron=1;
inputpars.ks_bandwidth=ks_ban;
inputpars.xlimtouse=[-0.5,5.5];
text(0,0.975*inputpars.yimtouse(2),['p grating r2 patt. comp. diff =',num2str(p_ttestr2_cpg)],'fontsize',12)
text(0,0.95*inputpars.yimtouse(2),['p plaid r2 patt. comp. diff =',num2str(p_ttestgr2_cpp)],'fontsize',12)
text(0,0.925*inputpars.yimtouse(2),['p gating plaid r2 diff comp. =',num2str(p_ttestpr2_comp)],'fontsize',12)
text(0,0.9*inputpars.yimtouse(2),['p gating plaid r2 diff patt. =',num2str(p_ttestpr2_patt)],'fontsize',12)
% plot violins
inputadata.inputdistrs=distribtouse;
inputpars.n_distribs=numel(inputadata.inputdistrs);
inputpars.dirstrcenters=(1:inputpars.n_distribs);
inputpars.xtickslabelvector{1}='grating - component';
inputpars.xtickslabelvector{2}='plaid - component';
inputpars.xtickslabelvector{3}='grating - pattern';
inputpars.xtickslabelvector{4}='plaid - pattern';
inputpars.distrcolors{1}=[50,200,0]./(1*255);
inputpars.distrcolors{2}=[50,200,0]./(2*255);
inputpars.distrcolors{3}=[255,150,0]./(1*255);
inputpars.distrcolors{4}=[255,150,0]./(2*255);
inputpars.boolscatteron=1;
inputaxh = plot_violinplot_PN_new(inputadata,inputpars); %#ok<NASGU>
xtickangle(45)
set(gca,'fontsize',12)
axis square
subplot(1,2,2)
hold on;
bar(1,nanmean(Explvar_T_comp_diff),'facecolor',inputpars.distrcolors{1},'edgecolor',inputpars.distrcolors{1})
bar(2,nanmean(Explvar_T_patt_diff),'facecolor',inputpars.distrcolors{3},'edgecolor',inputpars.distrcolors{3})
errorbar([1,2],...
    [nanmean(Explvar_T_comp_diff),nanmean(Explvar_T_patt_diff)],...
    [nanstd(Explvar_T_comp_diff)./sqrt(numel(Explvar_T_comp_diff)),nanstd(Explvar_T_patt_diff)./sqrt(numel(Explvar_T_patt_diff))],...
    'color','k','linewidth',2,'linestyle','none');
xlim([0,3])
ylimused=get(gca,'ylim');
text(0.25,0.95*ylimused(2),['p grating plaid patt. comp. diff = ',num2str(pp)],'fontsize',12);
ylabel('explained variance p vs. g difference')
xticks([1,2])
xticklabels({'component','pattern'})
title(['tuning explained variance p vs. g difference ( comp = ',...
    num2str(numel(distribtouse{1}),'%.0f'),...
    ' - patt  = ',num2str(numel(distribtouse{3}),'%.0f'),' )' ])
set(gca,'fontsize',12)
axis square
saveas(fighand1,[outfold,filesep,'RF_based_response_prediction_comparison_all_areas'],'jpg')
print(fighand1,'-depsc','-painters',[[outfold,filesep,'RF_based_response_prediction_comparison_all_areas'],'.eps'])

%% visualize tuning prediction results collapsed - r2

clear ylabellist yimtouselist ks_ban
distrtoplotlist{1}=collapsed_explvart_distribs;
distrtoplotlist{2}=collapsed_explvard_distribs;
ylabellist{1}='tuning explained variance';
ylabellist{2}='dynamics explained variance';
yimtouselist{1}=[-0.1,1];
yimtouselist{2}=[-0.1,1];
ks_ban{1}=0.1;
ks_ban{2}=0.1;
% initialize figure
fighand1=figure('units','normalized','outerposition',[0 0 1 1]);
for jj=1:2
    % decide wheter to use max or unrolled
    distribtouse=distrtoplotlist{jj}; % collapsed_max_fitted_rfs_r2_distribs
    subplot(1,2,jj)
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
    inputpars.titlestring=[ylabellist{jj},' ( P&G compRF = ',...
        num2str(numel(distribtouse{1}),'%.0f'),...
        ' - P&G pattRF  = ',num2str(numel(distribtouse{2}),'%.0f'),' )'];
    inputpars.boolscatteron=1;
    inputpars.ks_bandwidth=ks_ban{jj}; % inputpars.ks_bandwidth=0.25;
    %     inputpars.xlimtouse=[-0.5,5.5]; %[-1,5];
    inputpars.xlimtouse=[-0.5,4.5]; %[-1,5];
    % plot violins
    inputadata.inputdistrs=distribtouse;
    inputpars.n_distribs=numel(inputadata.inputdistrs);
    inputpars.dirstrcenters=(1:inputpars.n_distribs);
    %     inputpars.xtickslabelvector={'component - G','pattern - G','component - P','pattern - P'};
    inputpars.xtickslabelvector={'component','pattern','unclassified'};
    inputpars.distrcolors{1}=[50,200,0]./255;
    inputpars.distrcolors{2}=[255,150,0]./255;
    inputpars.distrcolors{3}=[0,0,0]./355;
    inputaxh = plot_violinplot_PN_new(inputadata,inputpars); %#ok<NASGU>
    pvalw = ranksum(distribtouse{1},distribtouse{2});
    [junk,pvalt] = ttest2(distribtouse{1},distribtouse{2});
    text(-0.25,0.95,['comp vs. pattern median diff p = ',num2str(pvalw)],'fontsize',12)
    set(gca,'fontsize',12)
    axis square
end
sgtitle(['RF- based response prediction analysis - all areas']) %#ok<NBRAK>
saveas(fighand1,[outfold,filesep,'RF_based_response_prediction_analysis_all_areas'],'jpg')
print(fighand1,'-depsc','-painters',[[outfold,filesep,'RF_based_response_prediction_analysis_all_areas'],'.eps'])

%% visualize tuning prediction results - partial correlation

% analyze prediction for consistent neurons only or not
bool_use_reclassified_labels=1;

% do preprocessing and select neurons to plot -----------------------------

% recompute observed cell type labels
pattcount=0;
recomputed_ctlabels=cell(1,3);
for cell_class_idx=1:length(recomputed_ctlabels)
    recomputed_ctlabels{cell_class_idx}=NaN(size(collapsed_Zc_O_distribs{cell_class_idx}));
    for i=1:numel(collapsed_Zc_O_distribs{cell_class_idx})
        if collapsed_Zp_O_distribs{cell_class_idx}(i)-max(collapsed_Zc_O_distribs{cell_class_idx}(i),0)>=1.28
            curr_ctlabel_O=2; % 2=pattern
            if cell_class_idx==2
                pattcount=pattcount+1;
            end
        elseif collapsed_Zc_O_distribs{cell_class_idx}(i)-max(collapsed_Zp_O_distribs{cell_class_idx}(i),0)>=1.28
            curr_ctlabel_O=1; % 1=component
        else
            curr_ctlabel_O=0; % 0=unclassified
        end
        recomputed_ctlabels{cell_class_idx}(i)=curr_ctlabel_O;
    end
end
% get consistent cells
idx_consist=cell(size(recomputed_ctlabels));
for cell_class_idx=1:length(recomputed_ctlabels)
    current_consist_idx=find(recomputed_ctlabels{cell_class_idx}==cell_types_codes(cell_class_idx));
    idx_consist{cell_class_idx}=current_consist_idx;
end
% compute predicted lables to compare with observed ones
recomputed_predicted_ctlabels=cell(1,3);
for cell_class_idx=1:length(recomputed_ctlabels)
    recomputed_predicted_ctlabels{cell_class_idx}=NaN(size(collapsed_Zc_P_distribs{cell_class_idx}));
    for i=1:numel(collapsed_Zc_P_distribs{cell_class_idx})
        if collapsed_Zp_P_distribs{cell_class_idx}(i)-max(collapsed_Zc_P_distribs{cell_class_idx}(i),0)>=1.28
            curr_ctlabel_O=2; % 2=pattern
            pattcount=pattcount+1;
        elseif collapsed_Zc_P_distribs{cell_class_idx}(i)-max(collapsed_Zp_P_distribs{cell_class_idx}(i),0)>=1.28
            curr_ctlabel_O=1; % 1=component
        else
            curr_ctlabel_O=0; % 0=unclassified
        end
        recomputed_predicted_ctlabels{cell_class_idx}(i)=curr_ctlabel_O;
    end
end
% restrict analysis to consistent cells - labels
recomputed_predicted_ctlabels_consist=cell(size(recomputed_predicted_ctlabels));
recomputed_ctlabels_consist=cell(size(recomputed_ctlabels));
for cell_class_idx=1:length(recomputed_predicted_ctlabels)
    current_consist_idx=idx_consist{cell_class_idx};
    recomputed_predicted_ctlabels_consist{cell_class_idx}=recomputed_predicted_ctlabels{cell_class_idx}(current_consist_idx);
    recomputed_ctlabels_consist{cell_class_idx}=recomputed_ctlabels{cell_class_idx}(current_consist_idx);
end
if not(bool_use_reclassified_labels)
    recomputed_ctlabels_to_use=recomputed_ctlabels;
    recomputed_ctlabels_predicted_to_use=recomputed_ctlabels;
else
    recomputed_ctlabels_to_use=recomputed_ctlabels_consist;
    recomputed_ctlabels_predicted_to_use=recomputed_predicted_ctlabels_consist;
end
% restrict analysis to consistent cells - Zc Zp obs and pred indexes
if not(bool_use_reclassified_labels)
    distrtoplotlistZcO=collapsed_Zc_O_distribs;
    distrtoplotlistZpO=collapsed_Zp_O_distribs;
    distrtoplotlistZcP=collapsed_Zc_P_distribs;
    distrtoplotlistZpP=collapsed_Zp_P_distribs;
else
    distrtoplotlistZcO=cell(1,length(recomputed_predicted_ctlabels));
    distrtoplotlistZpO=cell(1,length(recomputed_predicted_ctlabels));
    distrtoplotlistZcP=cell(1,length(recomputed_predicted_ctlabels));
    distrtoplotlistZpP=cell(1,length(recomputed_predicted_ctlabels));
    for cell_class_idx=1:length(recomputed_predicted_ctlabels)
        current_consist_idx=idx_consist{cell_class_idx};
        distrtoplotlistZcO{cell_class_idx}=collapsed_Zc_O_distribs{cell_class_idx}(current_consist_idx);
        distrtoplotlistZpO{cell_class_idx}=collapsed_Zp_O_distribs{cell_class_idx}(current_consist_idx);
        distrtoplotlistZcP{cell_class_idx}=collapsed_Zc_P_distribs{cell_class_idx}(current_consist_idx);
        distrtoplotlistZpP{cell_class_idx}=collapsed_Zp_P_distribs{cell_class_idx}(current_consist_idx);
    end
end

% plot Zp Zc observed predicted bullet  with PI comparison ----------------

fighand2=figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
axis square
hold on;
distrcolors{1}=[0.5,0.5,0.5];
distrcolors{2}=[50,200,0]./255;
distrcolors{3}=[255,150,0]./255;
for cell_class_idx=1:length(recomputed_ctlabels)
    for i=1:numel(distrtoplotlistZcO{cell_class_idx})
        curr_ctlabel_O=recomputed_ctlabels_to_use{cell_class_idx}(i);
        if curr_ctlabel_O~=0 && not(cell_class_idx==3)
            plot(distrtoplotlistZcO{cell_class_idx}(i),distrtoplotlistZpO{cell_class_idx}(i),'.','MarkerSize',45,'Color',distrcolors{curr_ctlabel_O+1});
            plot(distrtoplotlistZcP{cell_class_idx}(i),distrtoplotlistZpP{cell_class_idx}(i),'.','MarkerSize',45,'Color',distrcolors{curr_ctlabel_O+1}.*0.5);
            startvec=[distrtoplotlistZcO{cell_class_idx}(i),distrtoplotlistZpO{cell_class_idx}(i)];
            endvec=[distrtoplotlistZcP{cell_class_idx}(i),distrtoplotlistZpP{cell_class_idx}(i)];
            plot([startvec(1),endvec(1)],[startvec(2),endvec(2)],'linewidth',2,'Color',[distrcolors{curr_ctlabel_O+1}*0.75,0.25])
        end
    end
end
line([0 6], [1.28 7.28],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
line([1.28 7.28], [0 6],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
line([1.28 1.28], [-4 0],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
line([-4 0], [1.28 1.28],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
ylim([-4,7])
xlim([-4,7])
xlabel('Zc'); ylabel('Zp');
title(['partial correlation scatter (obs vs. pred) - all areas']); %#ok<NBRAK>
set(gca,'fontsize',12);
axis square
subplot(1,2,2)
axis square
hold on;
clear distrtoplotlist ylabellist yimtouselist ks_ban
% decide wheter to use full distribs or consistent enurons only for PI comparison
if not(bool_use_reclassified_labels)
    distrtoplotlist{1}=[collapsed_PI_O_distribs(1:2);collapsed_PI_P_distribs(1:2)];
elseif bool_use_reclassified_labels
    collapsed_PI_O_distribs_reclass=cell(1,length(recomputed_predicted_ctlabels));
    collapsed_PI_P_distribs_reclass=cell(1,length(recomputed_predicted_ctlabels));
    for cell_class_idx=1:length(recomputed_predicted_ctlabels)
        current_consist_idx=idx_consist{cell_class_idx};
        collapsed_PI_O_distribs_reclass{cell_class_idx}=collapsed_PI_O_distribs{cell_class_idx}(current_consist_idx);
        collapsed_PI_P_distribs_reclass{cell_class_idx}=collapsed_PI_P_distribs{cell_class_idx}(current_consist_idx);
    end
    distrtoplotlist{1}=[collapsed_PI_O_distribs_reclass(1:2)';collapsed_PI_P_distribs_reclass(1:2)'];
end
ylabellist{1}='pattern index';
yimtouselist{1}=[-7,7];
ks_ban{1}=0.65;
for jj=1
    % decide wheter to use max or unrolled
    distribtouse=distrtoplotlist{jj}; % collapsed_max_fitted_rfs_r2_distribs
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
    inputpars.titlestring=[ylabellist{jj},' ( comp # = ',...
        num2str(numel(distribtouse{1}),'%.0f'),...
        ' - patt #  = ',num2str(numel(distribtouse{2}),'%.0f'),' )'];
    inputpars.boolscatteron=1;
    inputpars.ks_bandwidth=ks_ban{jj}; % inputpars.ks_bandwidth=0.25;
    %     inputpars.xlimtouse=[-0.5,5.5]; %[-1,5];
    inputpars.xlimtouse=[-0.5,5.5]; %[-1,5];
    % plot violins
    inputadata.inputdistrs=distribtouse;
    inputpars.n_distribs=numel(inputadata.inputdistrs);
    inputpars.dirstrcenters=(1:inputpars.n_distribs);
    %     inputpars.xtickslabelvector={'component - G','pattern - G','component - P','pattern - P'};
    inputpars.xtickslabelvector={'component - obs ','pattern - obs','component - pred','pattern - pred'};
    inputpars.distrcolors{1}=[50,200,0]./255;
    inputpars.distrcolors{2}=[255,150,0]./255;
    inputpars.distrcolors{3}=[50,200,0]./355;
    inputpars.distrcolors{4}=[255,150,0]./355;
    inputaxh = plot_violinplot_PN_new(inputadata,inputpars); %#ok<NASGU>
    pvalw_comp = ranksum(distribtouse{1},distribtouse{3});
    pvalw_patt = ranksum(distribtouse{2},distribtouse{4});
    pvalw_displ = ranksum(distribtouse{1}-distribtouse{3},distribtouse{2}-distribtouse{4});
    text(-0.25,6.5,['comp median diff p = ',num2str(pvalw_comp)],'fontsize',12)
    text(-0.25,6.0,['patt median diff p = ',num2str(pvalw_patt)],'fontsize',12)
    text(-0.25,5.5,['patt vs. comp delta median diff p = ',num2str(pvalw_displ)],'fontsize',12)
    xtickangle(45)
    set(gca,'fontsize',12)
end
line([-0.5,5.5],[0,0],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
if bool_use_reclassified_labels
    sgtitle(['predicted vs. observed patternnes and componentness - all areas (consistent neurons only)']) %#ok<NBRAK>
else
    sgtitle(['predicted vs. observed patternnes and componentness - all areas']) %#ok<NBRAK>
end
saveas(fighand2,[outfold,filesep,'patt_comp_partial_correlation_obs_pred__all_areas'],'jpg')
print(fighand2,'-depsc','-painters',[[outfold,filesep,'patt_comp_partial_correlation_obs_pred__all_areas'],'.eps'])

% plot Zc Zp displacement comparison --------------------------------------

% set wheter to perform paired tests or not
pairedtestingbool=1;
% initialize figure
fighand2bis=figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
axis square
hold on;
clear distrtoplotlist ylabellist yimtouselist ks_ban
% decide wheter to use full distribs or consistent neurons only
if not(bool_use_reclassified_labels)
    distrtoplotlist{1}=[collapsed_Zc_O_distribs(1:2);collapsed_Zc_P_distribs(1:2)];
elseif bool_use_reclassified_labels
    collapsed_Zc_O_distribs_reclass=cell(1,length(recomputed_predicted_ctlabels));
    collapsed_Zc_P_distribs_reclass=cell(1,length(recomputed_predicted_ctlabels));
    for cell_class_idx=1:length(recomputed_predicted_ctlabels)
        current_consist_idx=idx_consist{cell_class_idx};
        collapsed_Zc_O_distribs_reclass{cell_class_idx}=collapsed_Zc_O_distribs{cell_class_idx}(current_consist_idx);
        collapsed_Zc_P_distribs_reclass{cell_class_idx}=collapsed_Zc_P_distribs{cell_class_idx}(current_consist_idx);
    end
    distrtoplotlist{1}=[collapsed_Zc_O_distribs_reclass(1:2)';collapsed_Zc_P_distribs_reclass(1:2)'];
end
ylabellist{1}='Zc';
yimtouselist{1}=[-4,8];
ks_ban{1}=0.65;
distrtoplotlist{1}=distrtoplotlist{1}([1,3,2,4]);
for jj=1
    % decide wheter to use max or unrolled
    distribtouse=distrtoplotlist{jj}; % collapsed_max_fitted_rfs_r2_distribs
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
    inputpars.titlestring=[ylabellist{jj},' ( comp # = ',...
        num2str(numel(distribtouse{1}),'%.0f'),...
        ' - patt #  = ',num2str(numel(distribtouse{2}),'%.0f'),' )'];
    inputpars.boolscatteron=1;
    inputpars.ks_bandwidth=ks_ban{jj}; % inputpars.ks_bandwidth=0.25;
    %     inputpars.xlimtouse=[-0.5,5.5]; %[-1,5];
    inputpars.xlimtouse=[-0.5,5.5]; %[-1,5];
    % plot violins
    inputadata.inputdistrs=distribtouse;
    inputpars.n_distribs=numel(inputadata.inputdistrs);
    inputpars.dirstrcenters=(1:inputpars.n_distribs);
    %     inputpars.xtickslabelvector={'component - G','pattern - G','component - P','pattern - P'};
    inputpars.xtickslabelvector={'component - obs ','component - pred','pattern - obs','pattern - pred'};
    inputpars.distrcolors{1}=[50,200,0]./255;
    inputpars.distrcolors{3}=[255,150,0]./255;
    inputpars.distrcolors{2}=[50,200,0]./355;
    inputpars.distrcolors{4}=[255,150,0]./355;
    [~,scatter_xs,scatter_ys] = plot_violinplot_PN_new(inputadata,inputpars);
    if not(pairedtestingbool)
        pvalw_comp = ranksum(distribtouse{1},distribtouse{2});
        pvalw_patt = ranksum(distribtouse{3},distribtouse{4});
        pvalw_displ = ranksum(distribtouse{1}-distribtouse{2},distribtouse{3}-distribtouse{4});
    else
        pvalw_comp = signrank(distribtouse{1},distribtouse{2});
        pvalw_patt = signrank(distribtouse{3},distribtouse{4});
        pvalw_displ = ranksum(distribtouse{1}-distribtouse{2},distribtouse{3}-distribtouse{4});
    end
    text(-0.25,6.5,['comp median diff p = ',num2str(pvalw_comp)],'fontsize',12)
    text(-0.25,6.0,['patt median diff p = ',num2str(pvalw_patt)],'fontsize',12)
    text(-0.25,5.5,['patt vs. comp delta median diff p = ',num2str(pvalw_displ)],'fontsize',12)
    xtickangle(45)
    set(gca,'fontsize',12)
end
hold on;
plot([scatter_xs{1},scatter_xs{2}]',[scatter_ys{1},scatter_ys{2}]','linewidth',2,'Color',[inputpars.distrcolors{1}*0.75,0.15])
plot([scatter_xs{3},scatter_xs{4}]',[scatter_ys{3},scatter_ys{4}]','linewidth',2,'Color',[inputpars.distrcolors{3}*0.75,0.15])
line([-0.5,5.5],[0,0],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
subplot(1,2,2)
axis square
hold on;
clear distrtoplotlist ylabellist yimtouselist ks_ban
% decide wheter to use full distribs or consistent neurons only
if not(bool_use_reclassified_labels)
    distrtoplotlist{1}=[collapsed_Zc_O_distribs(1:2);collapsed_Zp_P_distribs(1:2)];
elseif bool_use_reclassified_labels
    collapsed_Zp_O_distribs_reclass=cell(1,length(recomputed_predicted_ctlabels));
    collapsed_Zp_P_distribs_reclass=cell(1,length(recomputed_predicted_ctlabels));
    for cell_class_idx=1:length(recomputed_predicted_ctlabels)
        current_consist_idx=idx_consist{cell_class_idx};
        collapsed_Zp_O_distribs_reclass{cell_class_idx}=collapsed_Zp_O_distribs{cell_class_idx}(current_consist_idx);
        collapsed_Zp_P_distribs_reclass{cell_class_idx}=collapsed_Zp_P_distribs{cell_class_idx}(current_consist_idx);
    end
    distrtoplotlist{1}=[collapsed_Zp_O_distribs_reclass(1:2)';collapsed_Zp_P_distribs_reclass(1:2)'];
end
ylabellist{1}='Zp';
yimtouselist{1}=[-4,8];
distrtoplotlist{1}=distrtoplotlist{1}([1,3,2,4]);
ks_ban{1}=0.65;
for jj=1
    % decide wheter to use max or unrolled
    distribtouse=distrtoplotlist{jj}; % collapsed_max_fitted_rfs_r2_distribs
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
    inputpars.titlestring=[ylabellist{jj},' ( comp # = ',...
        num2str(numel(distribtouse{1}),'%.0f'),...
        ' - patt #  = ',num2str(numel(distribtouse{2}),'%.0f'),' )'];
    inputpars.boolscatteron=1;
    inputpars.ks_bandwidth=ks_ban{jj}; % inputpars.ks_bandwidth=0.25;
    %     inputpars.xlimtouse=[-0.5,5.5]; %[-1,5];
    inputpars.xlimtouse=[-0.5,5.5]; %[-1,5];
    % plot violins
    inputadata.inputdistrs=distribtouse;
    inputpars.n_distribs=numel(inputadata.inputdistrs);
    inputpars.dirstrcenters=(1:inputpars.n_distribs);
    %     inputpars.xtickslabelvector={'component - G','pattern - G','component - P','pattern - P'};
    inputpars.xtickslabelvector={'component - obs ','component - pred','pattern - obs','pattern - pred'};
    inputpars.distrcolors{1}=[50,200,0]./255;
    inputpars.distrcolors{3}=[255,150,0]./255;
    inputpars.distrcolors{2}=[50,200,0]./355;
    inputpars.distrcolors{4}=[255,150,0]./355;
    [~,scatter_xs,scatter_ys] = plot_violinplot_PN_new(inputadata,inputpars);
    if not(pairedtestingbool)
        pvalw_comp = ranksum(distribtouse{1},distribtouse{2});
        pvalw_patt = ranksum(distribtouse{3},distribtouse{4});
        pvalw_displ = ranksum(distribtouse{1}-distribtouse{2},distribtouse{3}-distribtouse{4});
        pvalw_compP_pattO = ranksum(distribtouse{2},distribtouse{3});
    else
        pvalw_comp = signrank(distribtouse{1},distribtouse{2});
        pvalw_patt = signrank(distribtouse{3},distribtouse{4});
        pvalw_displ = ranksum(distribtouse{1}-distribtouse{2},distribtouse{3}-distribtouse{4});
        pvalw_compP_pattO = ranksum(distribtouse{2},distribtouse{3});
    end
    text(-0.25,7,['compP pattO median diff p = ',num2str(pvalw_comp)],'fontsize',12)
    text(-0.25,6.5,['comp median diff p = ',num2str(pvalw_comp)],'fontsize',12)
    text(-0.25,6.0,['patt median diff p = ',num2str(pvalw_patt)],'fontsize',12)
    text(-0.25,5.5,['patt vs. comp delta median diff p = ',num2str(pvalw_displ)],'fontsize',12)
    xtickangle(45)
    set(gca,'fontsize',12)
end
hold on;
plot([scatter_xs{1},scatter_xs{2}]',[scatter_ys{1},scatter_ys{2}]','linewidth',2,'Color',[inputpars.distrcolors{1}*0.75,0.15])
plot([scatter_xs{3},scatter_xs{4}]',[scatter_ys{3},scatter_ys{4}]','linewidth',2,'Color',[inputpars.distrcolors{3}*0.75,0.15])
line([-0.5,5.5],[0,0],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
if bool_use_reclassified_labels
    sgtitle(['predicted vs. observed Zp and Zc - all areas (consistent neurons only)']) %#ok<NBRAK>
else
    sgtitle(['predicted vs. observed Zp and Zc - all areas']) %#ok<NBRAK>
end
saveas(fighand2bis,[outfold,filesep,'patt_comp_Zp_Zc_obs_pred_all_areas'],'jpg')
print(fighand2bis,'-depsc','-painters',[[outfold,filesep,'patt_comp_Zp_Zc_obs_pred_all_areas'],'.eps'])

%%
% barplot observed vs. predicted category change -------------------------.

fighand2tris=figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
reclass_ns_per_class=cell(1,2);
for cell_class_idx=1:2
    reclass_ns_per_class{cell_class_idx}=NaN(1,3);
    if cell_class_idx==1
        current_cell_type_code=1;
        current_cell_type_code_null=2;
    elseif cell_class_idx==2
        current_cell_type_code=2;
        current_cell_type_code_null=1;
    end
    reclass_ns_per_class{cell_class_idx}(1)=sum(recomputed_ctlabels_predicted_to_use{cell_class_idx}==current_cell_type_code);
    reclass_ns_per_class{cell_class_idx}(2)=sum(recomputed_ctlabels_predicted_to_use{cell_class_idx}==current_cell_type_code_null);
    reclass_ns_per_class{cell_class_idx}(3)=sum(recomputed_ctlabels_predicted_to_use{cell_class_idx}==0);
end
comp_n_dis=reclass_ns_per_class{1};
patt_n_dis=reclass_ns_per_class{2};
ytouse=[comp_n_dis./sum(comp_n_dis),patt_n_dis./sum(patt_n_dis)];
xtouse=[1:3,5:7];
hold on;
bar([xtouse([3]),0.1+xtouse([3])],[ytouse([3]),NaN],...
    'facecolor',[0.5,0.5,0.5],...
    'edgecolor',[0.5,0.5,0.5],...
    'facealpha',0.5,...
    'linewidth',1.5,...
    'barwidth',7)
bar([xtouse([6]),0.1+xtouse([6])],[ytouse([6]),NaN],...
    'facecolor',[0.5,0.5,0.5],...
    'edgecolor',[0.5,0.5,0.5],...
    'facealpha',0.5,...
    'linewidth',1.5,...
    'barwidth',7)
bar([xtouse([1]),0.1+xtouse([1])],[ytouse([1]),NaN],...
    'facecolor',compcol,...
    'edgecolor',compcol,...
    'facealpha',0.5,...
    'linewidth',1.5,...
    'barwidth',7)
bar([xtouse([4]),0.1+xtouse([4])],[ytouse([4]),NaN],...
    'facecolor',pattcol,...
    'edgecolor',pattcol,...
    'facealpha',0.5,...
    'linewidth',1.5,...
    'barwidth',7)
bar([xtouse([2]),0.1+xtouse([2])],[ytouse([2]),NaN],...
    'facecolor',compcol./2,...
    'edgecolor',compcol./2,...
    'facealpha',0.5,...
    'linewidth',1.5,...
    'barwidth',7)
bar([xtouse([5]),0.1+xtouse([5])],[ytouse([5]),NaN],...
    'facecolor',pattcol./2,...
    'edgecolor',pattcol./2,...
    'facealpha',0.5,...
    'linewidth',1.5,...
    'barwidth',7)
titlestring=['predicted vs. observed classification',' ( #comp = ',...
    num2str(sum(comp_n_dis),'%.0f'),...
    ' - #patt = ',num2str(sum(patt_n_dis),'%.0f'),' )'];
title(titlestring)
ylim([0,0.95])
xlim([0,8])
[chi2stat,p_chitest] = chiSquareTest([comp_n_dis;patt_n_dis]);
text(1,0.9,['chi square p = ',num2str(p_chitest)],'fontsize',12);
xticks(xtouse)
xticklabels({'same (cc)','opposite (cp)','unclassified (cu)','same (pp)','opposite (pc)','unclassified (pu)'})
xtickangle(45)
xlabel('')
ylabel('fraction of cells')
set(gca,'fontsize',12)
axis square
saveas(fighand2tris,[outfold,filesep,'pred_obs_classification_barplot_all'],'jpg')
print(fighand2tris,'-depsc','-painters',[[outfold,filesep,'pred_obs_classification_barplot_all'],'.eps'])

%% perform and visualize ibrahim analysis

% % % % set thresholds for ibrahim analysis
% % % r2_th_for_ibrh_anlaysis=0.4;
% % % dsi_th_for_ibrh_anlaysis=0; %0.2;
% % % dirs=0:30:330;
% % %
% % % % initialize data storage
% % % obs_grating_tc=cell(2,length(target_areas));
% % % obs_plaid_tc=cell(2,length(target_areas));
% % % pred_grating_tc=cell(2,length(target_areas));
% % % pred_plaid_tc=cell(2,length(target_areas));
% % % obs_grating_tc_realigned=cell(2,length(target_areas));
% % % pred_grating_tc_realigned=cell(2,length(target_areas));
% % % obs_plaid_tc_realigned=cell(2,length(target_areas));
% % % pred_plaid_tc_realigned=cell(2,length(target_areas));
% % % avg_filter_realigned=cell(2,length(target_areas));
% % %
% % % % loop over areas
% % % target_areas={'V1','LM','RL'};
% % % for target_area_idx=1:length(target_areas)
% % %
% % %     % set area
% % %     are=target_areas{target_area_idx};
% % %
% % %     % get nogocolor index by cell type and area
% % %     cell_types={'component','pattern'};
% % %     % loop over cell types
% % %     for cell_type_idx=1:2%length(cell_types)
% % %
% % %         % set cell type
% % %         cty=cell_types{cell_type_idx};
% % %
% % %         % load neuron selection results
% % %         sel_res_file=['D:\Backups\Personal_bk\PN_acute_analysis\processed_data',filesep,are,'_selected_neurons',filesep,'selected_population_',are,'.mat'];
% % %         load(sel_res_file)
% % %
% % %         % load fitted sta analysis results
% % %         sta_res_file=[outfold,filesep,'rf_fitting_results.mat'];
% % %         load(sta_res_file)
% % %
% % %         % select current category neurons with good enough Gabor fitted rf
% % %         current_r2s_per_frame=fitted_rfs_r2_distribs{cell_type_idx,target_area_idx};
% % %         orig_valid_neurons_idx=find(max(current_r2s_per_frame,[],1)>r2_th_for_pred_anlaysis);
% % %         current_r2s_per_frame=current_r2s_per_frame(:,orig_valid_neurons_idx);
% % %         valid_neurons_idx=find(max(current_r2s_per_frame,[],1)>r2_th_for_ibrh_anlaysis);
% % %         % TODO: fix DS selection
% % %         %         valid_ds_neu_idx=find(obs_DSI_distribs{cell_type_idx,target_area_idx}(valid_neurons_idx)>dsi_th_for_ibrh_anlaysis);
% % %         %         valid_neurons_idx=valid_neurons_idx(valid_ds_neu_idx);
% % %         current_selected_neu_idx=selected_cells_idx_distribs{cell_type_idx,target_area_idx}(valid_neurons_idx);
% % %         selectedsi_to_use=selectedsi(current_selected_neu_idx);
% % %         filters_to_use=fitted_rfs_distribs{cell_type_idx,target_area_idx}(:,:,:,valid_neurons_idx);
% % %         r2s_per_frame_to_use=current_r2s_per_frame(:,valid_neurons_idx);
% % %
% % %         % initialize tc datastructures
% % %         obs_grating_tc{cell_type_idx,target_area_idx}=nan(numel(current_selected_neu_idx),12);
% % %         obs_plaid_tc{cell_type_idx,target_area_idx}=nan(numel(current_selected_neu_idx),12);
% % %         pred_grating_tc{cell_type_idx,target_area_idx}=nan(numel(current_selected_neu_idx),12);
% % %         pred_plaid_tc{cell_type_idx,target_area_idx}=nan(numel(current_selected_neu_idx),12);
% % %         obs_grating_tc_realigned{cell_type_idx,target_area_idx}=nan(numel(current_selected_neu_idx),12);
% % %         pred_grating_tc_realigned{cell_type_idx,target_area_idx}=nan(numel(current_selected_neu_idx),12);
% % %         obs_plaid_tc_realigned{cell_type_idx,target_area_idx}=nan(numel(current_selected_neu_idx),12);
% % %         pred_plaid_tc_realigned{cell_type_idx,target_area_idx}=nan(numel(current_selected_neu_idx),12);
% % %         avg_filter_realigned{cell_type_idx,target_area_idx}=nan(size(filters_to_use,1),size(filters_to_use,2),numel(current_selected_neu_idx));
% % %
% % %         % loop over current good neurons
% % %         for i=1:numel(selectedsi_to_use)
% % %
% % %             %  get observed tuning curves
% % %             obs_grating_tc{cell_type_idx,target_area_idx}(i,:)=obs_COUNT_distribs{cell_type_idx,target_area_idx}(:,i,1);
% % %             obs_plaid_tc{cell_type_idx,target_area_idx}(i,:)=obs_COUNT_distribs{cell_type_idx,target_area_idx}(:,i,2);
% % %
% % %             % get observed average filter
% % %             current_goodframes=find(r2s_per_frame_to_use(:,i)<=r2_th_for_pred_anlaysis);
% % %             obs_avg_filter=nanmean(squeeze(filters_to_use(:,:,current_goodframes,i)),3);
% % %
% % %             % get temporary grating and plaid predicted tc
% % %             temp_obs_tc_grating=obs_grating_tc{cell_type_idx,target_area_idx}(i,:);
% % %             temp_obs_tc_plaid=obs_plaid_tc{cell_type_idx,target_area_idx}(i,:);
% % %             temp_pred_tc_grating=pred_COUNT_distribs{cell_type_idx,target_area_idx}(:,i,1);
% % %             temp_pred_tc_plaid=pred_COUNT_distribs{cell_type_idx,target_area_idx}(:,i,2);
% % %
% % %             % get predicted tuning curves
% % %             scaling_factor = max(max(temp_obs_tc_grating),max(temp_obs_tc_plaid))./...
% % %                 max(max(temp_pred_tc_grating),max(temp_pred_tc_plaid));
% % %             pred_grating_tc{cell_type_idx,target_area_idx}(i,:)=temp_pred_tc_grating.*scaling_factor;
% % %             pred_plaid_tc{cell_type_idx,target_area_idx}(i,:)=temp_pred_tc_plaid.*scaling_factor;
% % %
% % %             % reallign to tuning curve peak
% % %             [~,grating_idxmax]=randmax(obs_grating_tc{cell_type_idx,target_area_idx}(i,:));
% % %             obs_grating_tc_realigned{cell_type_idx,target_area_idx}(i,:)=circshift(obs_grating_tc{cell_type_idx,target_area_idx}(i,:),-(grating_idxmax-1));
% % %             pred_grating_tc_realigned{cell_type_idx,target_area_idx}(i,:)=circshift(pred_grating_tc{cell_type_idx,target_area_idx}(i,:),-(grating_idxmax-1));
% % %             obs_plaid_tc_realigned{cell_type_idx,target_area_idx}(i,:)=circshift(obs_plaid_tc{cell_type_idx,target_area_idx}(i,:),-(grating_idxmax-1));
% % %             pred_plaid_tc_realigned{cell_type_idx,target_area_idx}(i,:)=circshift(pred_plaid_tc{cell_type_idx,target_area_idx}(i,:),-(grating_idxmax-1));
% % %             avg_filter_realigned{cell_type_idx,target_area_idx}(:,:,i)=imresize(imrotate(obs_avg_filter,30*(grating_idxmax-1)),[size(obs_avg_filter,1),size(obs_avg_filter,2)]);
% % %
% % %             % TODO: plot average RF - rotation of image to work must have rf center in the center
% % %             %             figure;
% % %             %             subplot(1,3,1)
% % %             %             hold on;
% % %             %             plot(dirs,obs_grating_tc_realigned{cell_type_idx,target_area_idx}(i,:),'r')
% % %             %             plot(dirs,pred_grating_tc_realigned{cell_type_idx,target_area_idx}(i,:),'k')
% % %             %             subplot(1,3,3)
% % %             %             hold on;
% % %             %             plot(dirs,obs_grating_tc{cell_type_idx,target_area_idx}(i,:),'r')
% % %             %             plot(dirs,pred_grating_tc{cell_type_idx,target_area_idx}(i,:),'k')
% % %             %             subplot(1,3,2)
% % %             %             imagesc(obs_avg_filter)
% % %             %             title(num2str(current_selected_neu_idx(i)))
% % %         end
% % %
% % %     end
% % % end
% % % % initialize storage
% % % collapsed_obs_grating_tc=cell(2,1);
% % % collapsed_obs_plaid_tc=cell(2,1);
% % % collapsed_pred_grating_tc=cell(2,1);
% % % collapsed_pred_plaid_tc=cell(2,1);
% % % collapsed_obs_grating_tc_realigned=cell(2,1);
% % % collapsed_pred_grating_tc_realigned=cell(2,1);
% % % collapsed_obs_plaid_tc_realigned=cell(2,1);
% % % collapsed_pred_plaid_tc_realigned=cell(2,1);
% % % avg_obs_grating_tc_realigned=cell(2,1);
% % % avg_obs_plaid_tc_realigned=cell(2,1);
% % % avg_pred_grating_tc_realigned=cell(2,1);
% % % avg_pred_plaid_tc_realigned=cell(2,1);
% % % se_obs_grating_tc_realigned=cell(2,1);
% % % se_obs_plaid_tc_realigned=cell(2,1);
% % % se_pred_grating_tc_realigned=cell(2,1);
% % % se_pred_plaid_tc_realigned=cell(2,1);
% % % PI_avg_obs=cell(2,1);
% % % Zp_avg_obs=cell(2,1);
% % % Zc_avg_obs=cell(2,1);
% % % PI_avg_pred=cell(2,1);
% % % Zp_avg_pred=cell(2,1);
% % % Zc_avg_pred=cell(2,1);
% % % % loop over cell types
% % % for cell_type_idx=1:2
% % %     % loop over areas
% % %     for target_area_idx=1:length(target_areas)
% % %
% % %         collapsed_obs_grating_tc{cell_type_idx}=cat(1,collapsed_obs_grating_tc{cell_type_idx},obs_grating_tc{cell_type_idx,target_area_idx});
% % %         collapsed_obs_plaid_tc{cell_type_idx}=cat(1,collapsed_obs_plaid_tc{cell_type_idx},obs_plaid_tc{cell_type_idx,target_area_idx});
% % %         collapsed_pred_grating_tc{cell_type_idx}=cat(1,collapsed_pred_grating_tc{cell_type_idx},pred_grating_tc{cell_type_idx,target_area_idx});
% % %         collapsed_pred_plaid_tc{cell_type_idx}=cat(1,collapsed_pred_plaid_tc{cell_type_idx},pred_plaid_tc{cell_type_idx,target_area_idx});
% % %         collapsed_obs_grating_tc_realigned{cell_type_idx}=cat(1,collapsed_obs_grating_tc_realigned{cell_type_idx},obs_grating_tc_realigned{cell_type_idx,target_area_idx});
% % %         collapsed_obs_plaid_tc_realigned{cell_type_idx}=cat(1,collapsed_obs_plaid_tc_realigned{cell_type_idx},obs_plaid_tc_realigned{cell_type_idx,target_area_idx});
% % %         collapsed_pred_grating_tc_realigned{cell_type_idx}=cat(1,collapsed_pred_grating_tc_realigned{cell_type_idx},pred_grating_tc_realigned{cell_type_idx,target_area_idx});
% % %         collapsed_pred_plaid_tc_realigned{cell_type_idx}=cat(1,collapsed_pred_plaid_tc_realigned{cell_type_idx},pred_plaid_tc_realigned{cell_type_idx,target_area_idx});
% % %
% % %     end
% % % end
% % % % loop over cell types
% % % for cell_type_idx=1:2
% % %     avg_obs_grating_tc_realigned{cell_type_idx}=nanmean(collapsed_obs_grating_tc_realigned{cell_type_idx},1);
% % %     avg_obs_plaid_tc_realigned{cell_type_idx}=nanmean(collapsed_obs_plaid_tc_realigned{cell_type_idx},1);
% % %     avg_pred_grating_tc_realigned{cell_type_idx}=nanmean(collapsed_pred_grating_tc_realigned{cell_type_idx},1);
% % %     avg_pred_plaid_tc_realigned{cell_type_idx}=nanmean(collapsed_pred_plaid_tc_realigned{cell_type_idx},1);
% % %     se_obs_grating_tc_realigned{cell_type_idx}=nanstd(collapsed_obs_grating_tc_realigned{cell_type_idx},[],1)./sqrt(sum(not(isnan(collapsed_obs_grating_tc_realigned{cell_type_idx})),1));
% % %     se_obs_plaid_tc_realigned{cell_type_idx}=nanstd(collapsed_obs_plaid_tc_realigned{cell_type_idx},[],1)./sqrt(sum(not(isnan(collapsed_obs_plaid_tc_realigned{cell_type_idx})),1));
% % %     se_pred_grating_tc_realigned{cell_type_idx}=nanstd(collapsed_pred_grating_tc_realigned{cell_type_idx},[],1)./sqrt(sum(not(isnan(collapsed_pred_grating_tc_realigned{cell_type_idx})),1));
% % %     se_pred_plaid_tc_realigned{cell_type_idx}=nanstd(collapsed_pred_plaid_tc_realigned{cell_type_idx},[],1)./sqrt(sum(not(isnan(collapsed_pred_plaid_tc_realigned{cell_type_idx})),1));
% % %     [ PI_avg_obs{cell_type_idx}, ~, Zp_avg_obs{cell_type_idx}, Zc_avg_obs{cell_type_idx}, ~, ~, ~, ~, ~, ~ ]  = get_pattern_index( avg_obs_grating_tc_realigned{cell_type_idx}',avg_obs_plaid_tc_realigned{cell_type_idx}' );
% % %     [ PI_avg_pred{cell_type_idx}, ~, Zp_avg_pred{cell_type_idx}, Zc_avg_pred{cell_type_idx}, ~, ~, ~, ~, ~, ~ ]  = get_pattern_index( avg_pred_grating_tc_realigned{cell_type_idx}',avg_pred_plaid_tc_realigned{cell_type_idx}' );
% % % end
% % % % plot results
% % % fighand2=figure('units','normalized','outerposition',[0 0 1 1]);
% % % for cell_type_idx=1:2
% % %     subplot(1,2,cell_type_idx);
% % %     % plot polar plots with errorbars
% % %     hanldetouse=gca;
% % %     angvec=dirs;
% % %     alphatouse=0.10;
% % %     maxradiustouse=0.05*max([max(avg_obs_grating_tc_realigned{cell_type_idx}),max(avg_obs_plaid_tc_realigned{cell_type_idx})]);
% % %     linestyletouse='-';
% % %     muvec=avg_obs_grating_tc_realigned{cell_type_idx};
% % %     stdvec=se_obs_grating_tc_realigned{cell_type_idx};
% % %     colortouse=inputpars.distrcolors{cell_type_idx};
% % %     plot_polar_mu_std(hanldetouse,angvec,muvec,stdvec,colortouse,alphatouse,linestyletouse,maxradiustouse)
% % %     muvec=avg_obs_plaid_tc_realigned{cell_type_idx};
% % %     stdvec=se_obs_plaid_tc_realigned{cell_type_idx};
% % %     colortouse=inputpars.distrcolors{cell_type_idx}.*0.75;
% % %     plot_polar_mu_std(hanldetouse,angvec,muvec,stdvec,colortouse,alphatouse,linestyletouse,maxradiustouse)
% % %     linestyletouse='--';
% % %     muvec=avg_pred_grating_tc_realigned{cell_type_idx};
% % %     stdvec=se_pred_grating_tc_realigned{cell_type_idx};
% % %     colortouse=inputpars.distrcolors{cell_type_idx}.*0.5;
% % %     plot_polar_mu_std(hanldetouse,angvec,muvec,stdvec,colortouse,alphatouse,linestyletouse,maxradiustouse)
% % %     muvec=avg_pred_plaid_tc_realigned{cell_type_idx};
% % %     stdvec=se_pred_plaid_tc_realigned{cell_type_idx};
% % %     colortouse=inputpars.distrcolors{cell_type_idx}.*0.25;
% % %     plot_polar_mu_std(hanldetouse,angvec,muvec,stdvec,colortouse,alphatouse,linestyletouse,maxradiustouse)
% % %     ylimused=get(gca,'ylim');
% % %     ylimrange=diff(ylimused);
% % %     xlimused=get(gca,'xlim');
% % %     xlimrange=diff(xlimused);
% % %     text(hanldetouse,[xlimused(1)+0.05*xlimrange],0.9*[ylimused(1)+ylimrange],['obs PI = ',num2str(PI_avg_obs{cell_type_idx})],'fontsize',12) %#ok<NBRAK>
% % %     text(hanldetouse,[xlimused(1)+0.05*xlimrange],0.85*[ylimused(1)+ylimrange],['pred PI = ',num2str(PI_avg_pred{cell_type_idx})],'fontsize',12) %#ok<NBRAK>
% % %     text(hanldetouse,[xlimused(1)+0.05*xlimrange],0.75*[ylimused(1)+ylimrange],['pred Zp = ',num2str(Zp_avg_obs{cell_type_idx})],'fontsize',12) %#ok<NBRAK>
% % %     text(hanldetouse,[xlimused(1)+0.05*xlimrange],0.70*[ylimused(1)+ylimrange],['pred Zc = ',num2str(Zc_avg_pred{cell_type_idx})],'fontsize',12) %#ok<NBRAK>
% % %     set(gca,'fontsize',12)
% % %     title(['avg tc - ','all areas',' - ',cell_types{cell_type_idx},' ( # ',num2str(size(collapsed_obs_grating_tc_realigned{cell_type_idx},1)),' )'])
% % % end
% % % sgtitle(['realigned tuning curves analysis - all areas']) %#ok<NBRAK>
% % % saveas(fighand2,[outfold,filesep,'realigned_tuning_curves_analysis'],'jpg')

%% plot partial correlation plot - all areas

fighand4=figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
npatt_2=0;
ncomp_2=0;
inputpars.distrcolors{3}=[0,0,0];
inputpars.distrcolors{1}=[50,200,0]./255;
inputpars.distrcolors{2}=[255,150,0]./255;
% select indexes to use
switch labelstouse
    case 'oldbest'
        Zc_dist_touse=collapsed_Zc_O_oldbest_distribs;
        Zp_dist_touse=collapsed_Zp_O_oldbest_distribs;
    case 'old'
        Zc_dist_touse=collapsed_Zc_O_old_distribs;
        Zp_dist_touse=collapsed_Zp_O_old_distribs;
    case 'new'
        Zc_dist_touse=collapsed_Zc_O_distribs;
        Zp_dist_touse=collapsed_Zp_O_distribs;
    case 'oldbest_intersection'
        Zc_dist_touse=collapsed_Zc_O_oldbest_distribs;
        Zp_dist_touse=collapsed_Zp_O_oldbest_distribs;
end
% loop over cell classes
for cell_class_idx=1:numel(cell_types_codes)
    scatter(Zc_dist_touse{cell_class_idx},Zp_dist_touse{cell_class_idx},100,...
        'MarkerFaceColor',inputpars.distrcolors{cell_class_idx},'MarkerEdgeColor',inputpars.distrcolors{cell_class_idx},...
        'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
    valid_neu_idx=collapsed_obs_DSI_distribs{cell_class_idx}>=0.33;
    tempZp=Zp_dist_touse{cell_class_idx};
    tempZc=Zc_dist_touse{cell_class_idx};
    plot(tempZc(valid_neu_idx),tempZp(valid_neu_idx),'.','MarkerSize',25,'Color',inputpars.distrcolors{cell_class_idx});
    if cell_class_idx==1
        ncomp_2=sum(valid_neu_idx);
    elseif cell_class_idx==2
        npatt_2=sum(valid_neu_idx);
    else
        nuncl_2=sum(valid_neu_idx);
    end
end
ntot_2=npatt_2+ncomp_2+nuncl_2;
npatt_1=numel(Zc_dist_touse{2});
ncomp_1=numel(Zc_dist_touse{1});
nuncl_1=numel(Zc_dist_touse{3});
ntot_1=npatt_1+ncomp_1+nuncl_1;
line([0 6], [1.28 7.28],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
line([1.28 7.28], [0 6],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
line([1.28 1.28], [-4 0],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
line([-4 0], [1.28 1.28],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
xlabel('Zc'); ylabel('Zp');
title(['Partial correlation scatter - all areas',' ( ntot=',num2str(ntot_1),' ntotDS=',num2str(ntot_2),' )']);
text(-2.2,0.95*6.5,['n patt (DSI>0) = ',num2str(npatt_1),' ( frac = ',num2str(round(npatt_1/ntot_1,3)),' )'],'fontsize',12);
text(-2.2,0.91*6.5,['n comp (DSI>0) = ',num2str(ncomp_1),' ( frac = ',num2str(round(ncomp_1/ntot_1,3)),' )'],'fontsize',12);
text(-2.2,0.87*6.5,['n uncl (DSI>0) = ',num2str(nuncl_1),' ( frac = ',num2str(round(nuncl_1/ntot_1,3)),' )'],'fontsize',12);
text(-2.2,0.83*6.5,['n patt (DSI>0.33) = ',num2str(npatt_2),' ( frac = ',num2str(round(npatt_2/ntot_2,3)),' )'],'fontsize',12);
text(-2.2,0.79*6.5,['n comp (DSI>0.33) = ',num2str(ncomp_2),' ( frac = ',num2str(round(ncomp_2/ntot_2,3)),' )'],'fontsize',12);
text(-2.2,0.75*6.5,['n uncl (DSI>0.33) = ',num2str(nuncl_2),' ( frac = ',num2str(round(nuncl_2/ntot_2,3)),' )'],'fontsize',12);
set(gca,'fontsize',12);
axis square
ylim(gca,[-3.5,6.5])
xlim(gca,[-3.5,6.5])

axis square
saveas(fighand4,[outfold,filesep,'bullets_all_areas'],'jpg')
print(fighand4,'-depsc','-painters',[[outfold,filesep,'bullets_all_areas'],'.eps'])

%% plot partial correlation plot - each area separtely

areas={'V1','LM','RL'};
fighand4=figure('units','normalized','outerposition',[0 0 1 1]);
inputpars.distrcolors{3}=[0,0,0];
inputpars.distrcolors{1}=[50,200,0]./255;
inputpars.distrcolors{2}=[255,150,0]./255;
% select indexes to use
switch labelstouse
    case 'oldbest'
        Zc_dists_touse=Zc_O_oldbest_distribs;
        Zp_dists_touse=Zp_O_oldbest_distribs;
    case 'old'
        Zc_dists_touse=Zc_O_old_distribs;
        Zp_dists_touse=Zp_O_old_distribs;
    case 'new'
        Zc_dists_touse=Zc_O_distribs;
        Zp_dists_touse=Zp_O_distribs;
end
% loop over cell classes
for area_idx=1:3
    subplot(1,3,area_idx)
    hold on;
    npatt_2=0;
    ncomp_2=0;
    for cell_class_idx=1:numel(cell_types_codes)
        scatter(Zc_dists_touse{cell_class_idx,area_idx},Zp_dists_touse{cell_class_idx,area_idx},100,...
            'MarkerFaceColor',inputpars.distrcolors{cell_class_idx},'MarkerEdgeColor',inputpars.distrcolors{cell_class_idx},...
            'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
        valid_neu_idx=obs_DSI_distribs{cell_class_idx,area_idx}>=0.33;
        tempZp=Zp_dists_touse{cell_class_idx,area_idx};
        tempZc=Zc_dists_touse{cell_class_idx,area_idx};
        plot(tempZc(valid_neu_idx),tempZp(valid_neu_idx),'.','MarkerSize',25,'Color',inputpars.distrcolors{cell_class_idx});
        if cell_class_idx==1
            ncomp_2=sum(valid_neu_idx);
        elseif cell_class_idx==2
            npatt_2=sum(valid_neu_idx);
        else
            nuncl_2=sum(valid_neu_idx);
        end
    end
    ntot_2=npatt_2+ncomp_2+nuncl_2;
    npatt_1=numel(Zc_dists_touse{2,area_idx});
    ncomp_1=numel(Zc_dists_touse{1,area_idx});
    nuncl_1=numel(Zc_dists_touse{3,area_idx});
    ntot_1=npatt_1+ncomp_1+nuncl_1;
    line([0 6], [1.28 7.28],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
    line([1.28 7.28], [0 6],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
    line([1.28 1.28], [-4 0],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
    line([-4 0], [1.28 1.28],'LineWidth',1.5,'Color',[0.5,0.5,0.5]);
    xlabel('Zc'); ylabel('Zp');
    title(['Partial correlation scatter - ',areas{area_idx},' ( ntot=',num2str(ntot_1),' ntotDS=',num2str(ntot_2),' )']);
    text(-2.2,0.95*6.5,['n patt (DSI>0) = ',num2str(npatt_1),' ( frac = ',num2str(round(npatt_1/ntot_1,3)),' )'],'fontsize',12);
    text(-2.2,0.90*6.5,['n comp (DSI>0) = ',num2str(ncomp_1),' ( frac = ',num2str(round(ncomp_1/ntot_1,3)),' )'],'fontsize',12);
    text(-2.2,0.85*6.5,['n uncl (DSI>0) = ',num2str(nuncl_1),' ( frac = ',num2str(round(nuncl_1/ntot_1,3)),' )'],'fontsize',12);
    text(-2.2,0.80*6.5,['n patt (DSI>0.33) = ',num2str(npatt_2),' ( frac = ',num2str(round(npatt_2/ntot_2,3)),' )'],'fontsize',12);
    text(-2.2,0.75*6.5,['n comp (DSI>0.33) = ',num2str(ncomp_2),' ( frac = ',num2str(round(ncomp_2/ntot_2,3)),' )'],'fontsize',12);
    text(-2.2,0.70*6.5,['n uncl (DSI>0.33) = ',num2str(nuncl_2),' ( frac = ',num2str(round(nuncl_2/ntot_2,3)),' )'],'fontsize',12);
    set(gca,'fontsize',12);
    ylim(gca,[-3.5,6.5])
    xlim(gca,[-3.5,6.5])
    axis square
    comp_fracs1(area_idx)=ncomp_1/ntot_1; %#ok<SAGROW>
    patt_fracs1(area_idx)=npatt_1/ntot_1; %#ok<SAGROW>
    uncl_fracs1(area_idx)=nuncl_1/ntot_1; %#ok<SAGROW>
    comp_fracs2(area_idx)=ncomp_2/ntot_2; %#ok<SAGROW>
    patt_fracs2(area_idx)=npatt_2/ntot_2; %#ok<SAGROW>
    uncl_fracs2(area_idx)=nuncl_2/ntot_2; %#ok<SAGROW>
    ncomp_1s(area_idx)=ncomp_1; %#ok<SAGROW>
    npatt_1s(area_idx)=npatt_1; %#ok<SAGROW>
    nuncl_1s(area_idx)=nuncl_1; %#ok<SAGROW>
    ncomp_2s(area_idx)=ncomp_2; %#ok<SAGROW>
    npatt_2s(area_idx)=npatt_2; %#ok<SAGROW>
    nuncl_2s(area_idx)=nuncl_2; %#ok<SAGROW>
end
saveas(fighand4,[outfold,filesep,'bullets_all_areas_separately'],'jpg')
print(fighand4,'-depsc','-painters',[[outfold,filesep,'bullets_correlation_all_areas_separately'],'.eps'])

% compute ttest for difference of average PI per area
V1_Zc_distr=[Zc_dists_touse{1,1};Zc_dists_touse{2,1}];
LM_Zc_distr=[Zc_dists_touse{1,2};Zc_dists_touse{2,2}];
RL_Zc_distr=[Zc_dists_touse{1,3};Zc_dists_touse{2,3}];
V1_Zp_distr=[Zp_dists_touse{1,1};Zp_dists_touse{2,1}];
LM_Zp_distr=[Zp_dists_touse{1,2};Zp_dists_touse{2,2}];
RL_Zp_distr=[Zp_dists_touse{1,3};Zp_dists_touse{2,3}];
V1_PI_distr=V1_Zp_distr-V1_Zc_distr;
LM_PI_distr=LM_Zp_distr-LM_Zc_distr;
RL_PI_distr=RL_Zp_distr-RL_Zc_distr;
[p_V1_LM_PI,h_V1_LM_PI]=ttest2(V1_PI_distr,LM_PI_distr);
[p_V1_RL_PI,h_V1_RL_PI]=ttest2(V1_PI_distr,RL_PI_distr);
[p_RL_LM_PI,h_RL_LM_PI]=ttest2(RL_PI_distr,LM_PI_distr);

% initialize figure
fighand4bis = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
% get input for barplot
comp_fracs=comp_fracs1;
patt_fracs=patt_fracs1;
uncl_fracs=uncl_fracs1;
comp_xtouse=1:3;
patt_xtouse=1+4+(1:3);
uncl_xtouse=3+7+(1:3);
%get colors
compcol=[50,200,0]./255;
pattcol=[255,150,0]./255;
unclcol=[255/2,255/2,255/2]./255;
% plot bars
hold on;
bar(comp_xtouse,comp_fracs,...
    'facecolor',compcol,...
    'edgecolor',compcol,...
    'facealpha',0.5,...
    'linewidth',1.5...
    ) %#ok<*NBRAK>
bar(patt_xtouse,patt_fracs,...
    'facecolor',pattcol,...
    'edgecolor',pattcol,...
    'facealpha',0.5,...
    'linewidth',1.5...
    ) %#ok<*NBRAK>
bar(uncl_xtouse,uncl_fracs,...
    'facecolor',unclcol,...
    'edgecolor',unclcol,...
    'facealpha',0.5,...
    'linewidth',1.5...
    ) %#ok<*NBRAK>
titlestring=['all areas - fraction per class '];
title(titlestring)
ylim([0,0.7])
xlim([-1,uncl_xtouse(end)+2])
% [chi2stat,p_chitest] = chiSquareTest([comp_n_dis',patt_n_dis']);
% text(1,0.9,['chi square p = ',num2str(p_chitest)],'fontsize',12);
xtouse=[comp_xtouse,patt_xtouse,uncl_xtouse];
xticks(xtouse)
xtouselabls=[areas,areas,areas];
xticklabels(xtouselabls)
xtickangle(45)
xlabel('')
ylabel('fraction of cells')
set(gca,'fontsize',12)
axis square
subplot(1,2,2)
% get input for barplot
comp_fracs=comp_fracs2;
patt_fracs=patt_fracs2;
uncl_fracs=uncl_fracs2;
comp_xtouse=1:3;
patt_xtouse=1+4+(1:3);
uncl_xtouse=3+7+(1:3);
%get colors
compcol=0.75*[50,200,0]./255;
pattcol=0.75*[255,150,0]./255;
unclcol=0.75*[255/2,255/2,255/2]./255;
% plot bars
hold on;
bar(comp_xtouse,comp_fracs,...
    'facecolor',compcol,...
    'edgecolor',compcol,...
    'facealpha',0.5,...
    'linewidth',1.5...
    ) %#ok<*NBRAK>
bar(patt_xtouse,patt_fracs,...
    'facecolor',pattcol,...
    'edgecolor',pattcol,...
    'facealpha',0.5,...
    'linewidth',1.5...
    ) %#ok<*NBRAK>
bar(uncl_xtouse,uncl_fracs,...
    'facecolor',unclcol,...
    'edgecolor',unclcol,...
    'facealpha',0.5,...
    'linewidth',1.5...
    ) %#ok<*NBRAK>
titlestring=['all areas - fraction per class '];
title(titlestring)
ylim([0,0.7])
xlim([-1,uncl_xtouse(end)+2])
% [chi2stat,p_chitest] = chiSquareTest([comp_n_dis',patt_n_dis']);
% text(1,0.9,['chi square p = ',num2str(p_chitest)],'fontsize',12);
xtouse=[comp_xtouse,patt_xtouse,uncl_xtouse];
xticks(xtouse)
xtouselabls=[areas,areas,areas];
xticklabels(xtouselabls)
xtickangle(45)
xlabel('')
ylabel('fraction of cells')
set(gca,'fontsize',12)
axis square

% perform test for V1 - LM difference
tabl=[[npatt_2s(1),ncomp_2s(1)];[npatt_2s(2),ncomp_2s(2)]];
[chi2stat2,p_chitest2] = chiSquareTest(tabl);
tabl=[[npatt_1s(1),ncomp_1s(1)];[npatt_1s(2),ncomp_1s(2)]];
[chi2stat1,p_chitest1] = chiSquareTest(tabl);
x = table([npatt_1s(1:2)]',[ncomp_1s(1:2)]','VariableNames',{'p','c'},'RowNames',{'V1','LM'});
[~,p_frac_LM_V1_1,~] = fishertest(x,'Tail','both','Alpha',0.005);
x = table([npatt_2s(1:2)]',[ncomp_2s(1:2)]','VariableNames',{'p','c'},'RowNames',{'V1','LM'});
[~,p_frac_LM_V1_2,stats] = fishertest(x,'Tail','both','Alpha',0.005);

% save plot
saveas(fighand4bis,[outfold,filesep,'classification_barplots'],'jpg')
print(fighand4bis,'-depsc','-painters',[[outfold,filesep,'classification_barplots'],'.eps']);

%% re-plot single neuron examples

bool_replot_single_neurons_example=0;
bool_rerun_csi_analysis=0;
if bool_rerun_csi_analysis
    
    % load original tuning data
    tuning_file=['D:\Backups\Personal_bk\PN_acute_analysis\processed_data',filesep,'Tuning','.mat'];
    tuning=load(tuning_file);
    
    % get list of neurons to plot
    examples_all_n=[examples_comp_n_present;examples_patt_n_present];
    examples_all_ccode=[1*ones(size(examples_comp_n_present));2*ones(size(examples_patt_n_present))];
    examples_all_iidx=[1:length(examples_comp_n_present),1:length(examples_patt_n_present)]';
    
    % set list of conditions to plot for each neuron (pref SF and TF)
    ns_toplot=NaN(size(examples_all_n));
    SFs_toplot=NaN(size(examples_all_n));
    TFs_toplot=NaN(size(examples_all_n));
    for cell_class_idx=1:numel(cell_types_codes)
        for i=1:numel(examples_all_n)
            curr_selected_idx=find(collapsed_ns_distribs{cell_class_idx} == examples_all_n(i));
            if not(isnan(curr_selected_idx))
                SFs_toplot(i)=collapsed_pSF_distribs{cell_class_idx}(curr_selected_idx);
                TFs_toplot(i)=collapsed_pTF_distribs{cell_class_idx}(curr_selected_idx);
                ns_toplot(i)=examples_all_n(i);
            end
        end
    end
    
    % plot all single neuron examples
    for ii=1:numel(examples_all_n)
        % get current example neuron metadata
        i=find(examples_all_n==examples_all_n(ii));
        current_label=examples_all_ccode(i);
        input_sf=SFs_toplot(i);
        input_tf=TFs_toplot(i);
        input_neuron_num=ns_toplot(i);
        input_dirs=DIR;
        % set interp factor
        input_interp_factor=interp_factor;
        % get temporal parameters
        sr=videoframerate;
        countingwindowlim=[0,1];
        % get observed psths
        obs_FR=NaN(sr*diff(countingwindowlim),length(input_dirs),2);
        obs_COUNT=NaN(1,length(input_dirs),2);
        obs_RASTER=cell(1,length(input_dirs),2);
        pDIR=cell(1,2);
        for k=1:2
            input_stimtype=stimulustypes{k};
            for dir_idx=1:length(input_dirs)
                tic
                % get current direction
                current_direction=input_dirs(dir_idx);
                % get observed psth
                [current_psth,input_RASTER,current_psth_edges] = get_psth_PN(input_neuron_num,input_sf,input_tf,current_direction,input_stimtype);
                % get counting window span
                counting_window_span=countingwindowlim(2)-countingwindowlim(1);
                % smooth observed psth to match timescale of prediction
                gaussFilter = gausswin(10,1/0.5);
                gaussFilter = gaussFilter / sum(gaussFilter);
                current_psth=conv(current_psth, gaussFilter,'same');
                % resample observed psth
                current_psth=squeeze(current_psth);
                new_samplenum=sr.*(counting_window_span);
                old_samplenum=length(current_psth);
                current_psth_resampled = resample(current_psth,new_samplenum,old_samplenum);
                % get observed variables
                obs_FR(:,dir_idx,k)=current_psth_resampled;
                obs_COUNT(dir_idx,k)=nansum(current_psth_resampled);
                obs_RASTER{dir_idx,k}=input_RASTER;
                toc
            end
            % get preferred direction
            pDIR{k}=DIR(find(obs_COUNT(:,k)==max(obs_COUNT(:,k))));
        end
        % recover original Zp and Zc
        i_collapsed_id=examples_idx_present{examples_all_ccode(i)}(examples_all_iidx(i));
        original_dsi=collapsed_dsi_distribs{examples_all_ccode(i)}(i_collapsed_id);
        original_Zc=collapsed_Zc_O_oldbest_distribs{examples_all_ccode(i)}(i_collapsed_id);
        original_Zp=collapsed_Zp_O_oldbest_distribs{examples_all_ccode(i)}(i_collapsed_id);
        % recompute Zp and Zc
        tuning_curve_grating=obs_COUNT(:,1);
        tuning_curve_plaid=obs_COUNT(:,2);
        [ ~, ~, recomputed_Zp, recomputed_Zc, ~, ~, ~, ~, ~, ~ ] =...
            get_pattern_index( tuning_curve_grating,tuning_curve_plaid );
        % to check from scratch originally computed Zp and Zc values
        %     tuning.Zc(ns_toplot(i),SFs_toplot(i)==SF,TFs_toplot(i)==TF)
        %     tuning.Zp(ns_toplot(i),SFs_toplot(i)==SF,TFs_toplot(i)==TF)
        %     tcg_to_plot=tuning.tuning_curve_z(:,ns_toplot(i),SFs_toplot(i)==SF,TFs_toplot(i)==TF,1);
        %     tcp_to_plot=tuning.tuning_curve_z(:,ns_toplot(i),SFs_toplot(i)==SF,TFs_toplot(i)==TF,2);
        %     [ ~, ~, temp1, temp2, ~, ~, ~, ~, ~, ~ ] =...
        %         get_pattern_index( tcg_to_plot,tcp_to_plot );
        % plot figure ----------
        ffffffff = figure('units','normalized','outerposition',[0 0 1 1]);
        % set psth subplot position
        sb1=subplot(1,2,2);
        axis square
        % get psths to plot
        pipi=NaN(1,2);
        for k=1:2
            hold on;
            psth_observed=obs_FR(:,DIR==pDIR{1},k);
            psth_observed=psth_observed./(max(psth_observed(:)));
            % set color and tag to use
            if current_label==2
                coltuse=[255,150,0]./255;
            elseif current_label==1
                coltuse=[50,200,0]./255;
            elseif current_label==0
                coltuse=[150,150,150]./255;
            end
            hold on;
            % draw psth
            pipi(k)=plot(gca,(0:length(psth_observed)).*(1/sr),[psth_observed(1);psth_observed],'-','Color',coltuse./(k),'LineWidth',2.5);
            plot_shaded_auc(gca,(0:length(psth_observed))*(1/sr),[psth_observed(1);psth_observed]',0.15,coltuse./(k))
        end
        % draw raster
        for k=1:2
            raster=obs_RASTER{pDIR{1}==DIR,k};
            spp=raster;
            spk=0;
            for kkk=1:numel(raster)
                if not(isempty(raster{kkk}))
                    plot(raster{kkk},1.05+0.05*kkk,'.','color',coltuse./(k),'MarkerSize',15)
                    spk=spk+numel(raster{kkk});
                else
                end
            end
            xlim([-0.2,1.2]);
            ylim([-0,4]);
            if k==1
                tt=text(0.05,3.8,['pDIR grating = ',num2str(pDIR{1}),' d'],'FontSize',12);
                ttt=text(0.05,3.6,['spike count grating = ',num2str(spk)],'FontSize',12);
            else
                tt2=text(0.05,3.4,['pDIR plaid = ',num2str(pDIR{2}),' d'],'FontSize',12);
                ttt2=text(0.05,3.2,['spike count plaid = ',num2str(spk)],'FontSize',12);
            end
        end
        plot([0,0],[0,5],'--k', 'LineWidth',2)
        plot([1,1],[0,5],'--k', 'LineWidth',2)
        hlabelx=get(gca,'Xlabel');
        set(hlabelx,'String','time (s)','FontSize',12,'color','k')
        hlabely=get(gca,'Ylabel');
        set(hlabely,'String','normalized firing rate','FontSize',12,'color','k')
        legend(gca,pipi,{'grating','plaid'})
        title(['raster and psth (n=',num2str(input_neuron_num),') - DIR=',num2str(pDIR{1}),' - TF=',num2str(input_tf),' - SF=',num2str(input_sf)])
        set(gca,'FontSize',12);
        % set psth subplot position
        ppol=polaraxes('Position',[-0.05,0.19,.65,.65]);
        hold on;
        for k=1:2
            temp_obs_tc=squeeze(obs_COUNT(:,k))'./max(squeeze(obs_COUNT(:,k)));
            obs_tc=[temp_obs_tc,temp_obs_tc(1)];
            % draw plar plots
            p2=polarplot(ppol,[deg2rad(DIR),2*pi],obs_tc,'-');
            set(p2,'color',coltuse./k)
            set(p2, 'linewidth', 3.5);
        end
        title(ppol,['polar plots (n=',num2str(input_neuron_num),') - DIR=',num2str(pDIR{1}),' - TF=',num2str(input_tf),' - SF=',num2str(input_sf)])
        tx2=text(ppol,deg2rad(45),1.15,['orig DSI = ',num2str(original_dsi,'%.01f')],'fontsize',12);
        tx3=text(ppol,deg2rad(40),1.15,['original Zc = ',num2str(original_Zc,'%.01f')],'fontsize',12);
        tx4=text(ppol,deg2rad(35),1.15,['recomputed Zc = ',num2str(recomputed_Zc,'%.01f')],'fontsize',12);
        tx5=text(ppol,deg2rad(30),1.15,['original Zp = ',num2str(original_Zp,'%.01f')],'fontsize',12);
        tx6=text(ppol,deg2rad(25),1.15,['recomputed Zp = ',num2str(recomputed_Zp,'%.01f')],'fontsize',12);
        set(ppol,'fontsize',12);
        % save figure
        saveas(ffffffff,[outfold,filesep,'single_neu_example_pSF_pTF_n',num2str(input_neuron_num)],'jpg')
        print(ffffffff,'-depsc','-painters',[[outfold,filesep,'single_neu_example_pSF_pTF_n',num2str(input_neuron_num)],'.eps'])
        close all
    end
end

%% analyze plaid vs. grating responses at their best ----------------------

bool_rerun_csi_analysis=0;
bool_plot_csi_diagnostics=0;
if not(bool_rerun_csi_analysis) && exist(['CSI_results_',labelstouse,'.mat']) %#ok<EXIST>
    
    % load previously computaed CSI
    load(['CSI_results_',labelstouse,'.mat'])
    
elseif bool_rerun_csi_analysis
    
    % load original tuning data
    tuning_file=['D:\Backups\Personal_bk\PN_acute_analysis\processed_data',filesep,'Tuning','.mat']; %#ok<*UNRCH>
    tuning=load(tuning_file);
    
    % get parameters
    sr=videoframerate;
    countingwindowlim=[0,1];
    input_dirs=DIR;
    
    % set interp factor
    input_interp_factor=interp_factor;
    
    % set best spike count and CSI structures
    best_spike_counts{cell_class_idx}=cell(1,numel(cell_types_codes));
    CSI{cell_class_idx}=cell(1,numel(cell_types_codes));
    
    % loop over cell classes
    for cell_class_idx=1:numel(cell_types_codes)
        
        % get neurons
        current_ns=collapsed_ns_distribs{cell_class_idx};
        current_g_pSFs=pref_SF(current_ns,1);
        current_g_pTFs=pref_TF(current_ns,1);
        current_p_pSFs=pref_SF(current_ns,1);
        current_p_pTFs=pref_TF(current_ns,1);
        
        % initialize best spike count and CSI matrices
        best_spike_counts{cell_class_idx}=NaN(2,numel(current_ns));
        CSI{cell_class_idx}=NaN(1,numel(current_ns));
        
        % loop over neurons
        for i=1:numel(current_ns)
            
            % initalize observed psthsh structures
            obs_FR_adapted=NaN(sr*diff(countingwindowlim),length(input_dirs),numel(stimulustypes));
            obs_COUNT_adapted=NaN(1,length(input_dirs),numel(stimulustypes));
            obs_RASTER_adapted=cell(1,length(input_dirs),numel(stimulustypes));
            pDIR_adapted=cell(1,numel(stimulustypes));
            % recover psths of best condition and pref dir
            for k=1:numel(stimulustypes)
                current_label=cell_class_idx;
                if k==1
                    input_sf=current_g_pSFs(i);
                    input_tf=current_g_pTFs(i);
                elseif k==2
                    input_sf=current_p_pSFs(i);
                    input_tf=current_p_pTFs(i);
                end
                input_neuron_num=current_ns(i);
                input_dirs=DIR;
                input_stimtype=stimulustypes{k};
                for dir_idx=1:length(input_dirs)
                    tic
                    % get current direction
                    current_direction=input_dirs(dir_idx);
                    % get observed psth
                    [current_psth,input_RASTER,current_psth_edges] = get_psth_PN(input_neuron_num,input_sf,input_tf,current_direction,input_stimtype);
                    % get counting window span
                    counting_window_span=countingwindowlim(2)-countingwindowlim(1);
                    % smooth observed psth to match timescale of prediction
                    gaussFilter = gausswin(10,1/0.5);
                    gaussFilter = gaussFilter / sum(gaussFilter);
                    current_psth=conv(current_psth, gaussFilter,'same');
                    % resample observed psth
                    current_psth=squeeze(current_psth);
                    new_samplenum=sr.*(counting_window_span);
                    old_samplenum=length(current_psth);
                    current_psth_resampled = resample(current_psth,new_samplenum,old_samplenum);
                    % get observed variables
                    obs_FR(:,dir_idx,k)=current_psth_resampled;
                    obs_COUNT(dir_idx,k)=nansum(current_psth_resampled);
                    obs_RASTER{dir_idx,k}=input_RASTER;
                    toc
                end
                % get preferred direction
                pDIR{k}=DIR(find(obs_COUNT(:,k)==max(obs_COUNT(:,k))));
            end
            % compute best spike counts
            for k=1:2
                raster=obs_RASTER{pDIR{k}==DIR,k};
                spp=raster;
                spk=0;
                for kkk=1:numel(raster)
                    if not(isempty(raster{kkk}))
                        spk=spk+numel(raster{kkk});
                    else
                    end
                end
                best_spike_counts{cell_class_idx}(k,i)=spk;
            end
            % compute CSI as normalized difference
            CSI{cell_class_idx}(i)=(best_spike_counts{cell_class_idx}(1,i)-...
                best_spike_counts{cell_class_idx}(2,i))./...
                nansum(best_spike_counts{cell_class_idx}(:,i));
            
            if bool_plot_csi_diagnostics
                % plot figure
                fffffff = figure('units','normalized','outerposition',[0 0 1 1]);
                % set psth subplot position
                sb1=subplot(1,2,2);
                axis square
                % get psths to plot
                pipi=NaN(1,2);
                for k=1:2
                    hold on;
                    psth_observed=obs_FR(:,DIR==pDIR{1},k);
                    psth_observed=psth_observed./(max(psth_observed(:)));
                    % set color and tag to use
                    if current_label==2
                        coltuse=[255,150,0]./255;
                    elseif current_label==1
                        coltuse=[50,200,0]./255;
                    elseif current_label==0
                        coltuse=[150,150,150]./255;
                    end
                    hold on;
                    % draw psth
                    pipi(k)=plot(gca,(0:length(psth_observed)).*(1/sr),[psth_observed(1);psth_observed],'-','Color',coltuse./(k),'LineWidth',2.5);
                    plot_shaded_auc(gca,(0:length(psth_observed))*(1/sr),[psth_observed(1);psth_observed]',0.15,coltuse./(k))
                end
                % draw raster
                for k=1:2
                    raster=obs_RASTER{pDIR{k}==DIR,k};
                    spp=raster;
                    spk=0;
                    for kkk=1:numel(raster)
                        if not(isempty(raster{kkk}))
                            plot(raster{kkk},1.05+0.05*kkk,'.','color',coltuse./(k),'MarkerSize',15)
                            spk=spk+numel(raster{kkk});
                        else
                        end
                    end
                    xlim([-0.2,1.2]);
                    ylim([-0,4]);
                    if k==1
                        tt=text(0.05,3.8,['pDIR grating = ',num2str(pDIR{1}),' d'],'FontSize',12);
                        ttt=text(0.05,3.6,['spike count grating = ',num2str(spk)],'FontSize',12);
                    else
                        tt2=text(0.05,3.4,['pDIR plaid = ',num2str(pDIR{2}),' d'],'FontSize',12);
                        ttt2=text(0.05,3.2,['spike count plaid = ',num2str(spk)],'FontSize',12);
                    end
                end
                plot([0,0],[0,5],'--k', 'LineWidth',2)
                plot([1,1],[0,5],'--k', 'LineWidth',2)
                hlabelx=get(gca,'Xlabel');
                set(hlabelx,'String','time (s)','FontSize',12,'color','k')
                hlabely=get(gca,'Ylabel');
                set(hlabely,'String','normalized firing rate','FontSize',12,'color','k')
                legend(gca,pipi,{'grating','plaid'})
                title(['raster and psth (n=',num2str(input_neuron_num),') - DIR=','best',' - TF=',num2str(input_tf),' - SF=',num2str(input_sf)])
                set(gca,'FontSize',12);
                % set psth subplot position
                ppol=polaraxes('Position',[-0.05,0.19,.65,.65]);
                hold on;
                for k=1:2
                    temp_obs_tc=squeeze(obs_COUNT(:,k))'./max(squeeze(obs_COUNT(:,k)));
                    obs_tc=[temp_obs_tc,temp_obs_tc(1)];
                    % draw plar plots
                    p2=polarplot(ppol,[deg2rad(DIR),2*pi],obs_tc,'-');
                    set(p2,'color',coltuse./k)
                    set(p2, 'linewidth', 3.5);
                    % draw polarscatter of pref dir
                    ps2=polarscatter(ppol,deg2rad(pDIR{k}),obs_tc(find(DIR==pDIR{k})),...
                        125,'markerfacecolor',coltuse./k,'markeredgecolor',coltuse./k);
                end
                title(ppol,['polar plots (n=',num2str(input_neuron_num),') - DIR=',num2str(pDIR{1}),' - TF=',num2str(input_tf),' - SF=',num2str(input_sf)])
                tx2=text(ppol,deg2rad(45),1.15,['CSI = ',num2str(CSI{cell_class_idx}(i),'%.01f')],'fontsize',12);
                set(ppol,'fontsize',12);
                axis square
                % save figure
                saveas(fffffff,[outfold,filesep,'CSI_diagnostics_n',num2str(input_neuron_num)],'jpg')
                print(fffffff,'-depsc','-painters',[[outfold,filesep,'CSI_diagnostics_n',num2str(input_neuron_num)],'.eps'])
                close all
            end
            
            fprintf(['neuron#',num2str(i),' class#',num2str(cell_class_idx),' analyzed...\n'])
            
        end
    end
    save(['CSI_results_',labelstouse,'.mat'])
end

% get CSI distributions restricted to DS units only
CSI_ds=cell(1,numel(CSI));
OSI_ds=cell(1,numel(CSI));
DSI_ds=cell(1,numel(CSI));
ns_ds=cell(1,numel(CSI));
for cell_type_idx=1:numel(CSI)
    selected_ds_idxs=collapsed_obs_DSI_distribs{cell_type_idx}>=0.33;
    CSI_ds{cell_type_idx}=CSI{cell_type_idx}(selected_ds_idxs);
    OSI_ds{cell_type_idx}=collapsed_obs_OSI_distribs{cell_type_idx}(selected_ds_idxs);
    DSI_ds{cell_type_idx}=collapsed_obs_DSI_distribs{cell_type_idx}(selected_ds_idxs);
    ns_ds{cell_type_idx}=collapsed_ns_distribs{cell_type_idx}(selected_ds_idxs);
end

% decide wheter to plot full or DS distribution
boolplotCSIofDSonly=0;
if boolplotCSIofDSonly
    CSI_to_use=CSI_ds;
    nametag='CSIds';
    DSI_to_use=DSI_ds;
    nametagbis='DSIds';
    ns_to_use=ns_ds;
else
    CSI_to_use=CSI;
    nametag='CSIds';
    DSI_to_use=collapsed_obs_DSI_distribs;
    nametagbis='DSIds';
    ns_to_use=collapsed_ns_distribs;
end

% load original indexing data to get CSI per area
indexing_file=['D:\Backups\Personal_bk\PN_acute_analysis\processed_data',filesep,'Indexing','.mat']; %#ok<*UNRCH>
indexing=load(indexing_file);
% initialize CSI per area
CSI_to_use_per_area=cell(numel(areas),1);
area_codes=[0,1,2];
area_labels={'V1','LM','RL'};
% loop over areas
for area_idx=1:numel(areas)
    % get current areacode
    selected_areacode=area_codes(area_idx);
    % initialize CSI per area structure
    CSI_to_use_per_area{area_idx}=[];
    % loop over cell types
    for cell_type_idx=1:numel(CSI)
        % get current cell type areacodes
        selected_ns=ns_to_use{cell_type_idx};
        selected_arealabels=M(selected_ns,4);
        % select current neuron idxes
        current_neu_idx=find(selected_arealabels==selected_areacode);
        if not(isempty(current_neu_idx))
            CSI_to_use_per_area{area_idx}=[CSI_to_use_per_area{area_idx},...
                CSI_to_use{cell_type_idx}(current_neu_idx)];
        end
    end
end


% plot CSI analysis result
fighand0=figure('units','normalized','outerposition',[0 0 1 1]);
% plot distributions of CSI ---
subplot(1,2,1)
% decide wheter to use max or unrolled
distribtouse=CSI_to_use(1:3); % collapsed_max_fitted_rfs_r2_distribs
inputpars.inputaxh=gca;
hold(inputpars.inputaxh,'on')
% set settings for violin distribution plotting
inputpars.boxplotwidth=0.4;%0.5;
inputpars.boxplotlinewidth=2;
inputpars.densityplotwidth=0.5;%0.5;
inputpars.yimtouse=[-0.75,0.75];
% inputpars.yimtouse=[0,8];
inputpars.scatterjitter=inputpars.boxplotlinewidth*0.1;
inputpars.scatteralpha=0.15;
inputpars.scattersize=20;
inputpars.distralpha=0.5;
inputpars.xlabelstring=[];
inputpars.ylabelstring='CSI';
inputpars.titlestring=['cross-suppression index (',nametag,')',' ( comp # = ',...
    num2str(numel(distribtouse{1}),'%.0f'),...
    ' - patt #  = ',num2str(numel(distribtouse{2}),'%.0f'),...
    ' - uncl #  = ',num2str(numel(distribtouse{3}),'%.0f'),' )'];
inputpars.boolscatteron=1;
inputpars.ks_bandwidth=0.05;
inputpars.xlimtouse=[-0.5,4.5]; %[-1,5];
% plot violins
inputadata.inputdistrs=distribtouse;
inputpars.n_distribs=numel(inputadata.inputdistrs);
inputpars.dirstrcenters=(1:inputpars.n_distribs);
inputpars.xtickslabelvector={'component','pattern','unclassified'};
inputpars.distrcolors{1}=[50,200,0]./255;
inputpars.distrcolors{2}=[255,150,0]./255;
inputaxh = plot_violinplot_PN_new(inputadata,inputpars); %#ok<NASGU>
% [~,pval_csi_pc] = ttest2(distribtouse{1}',distribtouse{2}');
% [~,pval_csi_c] = ttest(distribtouse{1});
% [~,pval_csi_p] = ttest(distribtouse{2});
% [~,pval_csi_unc] = ttest(distribtouse{3});
[pval_csi_pc,~] = ranksum(distribtouse{1}',distribtouse{2}');
[pval_csi_c,~] = signrank(distribtouse{1});
[pval_csi_p,~] = signrank(distribtouse{2});
[pval_csi_unc,~] = signrank(distribtouse{3});
usedxlim=get(gca,'xlim');
hold on;
plot(gca,usedxlim,[0,0],'--','linewidth',2,'color',[0.5,0.5,0.5])
text(gca,-0.25,-0.5,['comp vs. patt median csi diff p = ',num2str(pval_csi_pc)],'fontsize',12)
text(gca,-0.25,-0.55,['comp median csi p = ',num2str(pval_csi_c)],'fontsize',12)
text(gca,-0.25,-0.6,['patt median csi p = ',num2str(pval_csi_p)],'fontsize',12)
text(gca,-0.25,-0.65,['uncl median csi p = ',num2str(pval_csi_unc)],'fontsize',12)
xtickangle(45)
set(gca,'fontsize',12)
axis square
% % plot distributions of DSI ---
% subplot(1,2,2)
% % decide wheter to use max or unrolled
% distribtouse=DSI_to_use(1:3); % collapsed_max_fitted_rfs_r2_distribs
% inputpars.inputaxh=gca;
% hold(inputpars.inputaxh,'on')
% % set settings for violin distribution plotting
% inputpars.boxplotwidth=0.4;%0.5;
% inputpars.boxplotlinewidth=2;
% inputpars.densityplotwidth=0.5;%0.5;
% inputpars.yimtouse=[-0.75,1.05];
% % inputpars.yimtouse=[0,8];
% inputpars.scatterjitter=inputpars.boxplotlinewidth*0.1;
% inputpars.scatteralpha=0.15;
% inputpars.scattersize=20;
% inputpars.distralpha=0.5;
% inputpars.xlabelstring=[];
% inputpars.ylabelstring='DSI';
% inputpars.titlestring=['dsi (',nametagbis,')',' ( comp # = ',...
%     num2str(numel(distribtouse{1}),'%.0f'),...
%     ' - patt #  = ',num2str(numel(distribtouse{2}),'%.0f'),...
%     ' - uncl #  = ',num2str(numel(distribtouse{3}),'%.0f'),' )'];
% inputpars.boolscatteron=1;
% inputpars.ks_bandwidth=0.05;
% inputpars.xlimtouse=[-0.5,4.5]; %[-1,5];
% % plot violins
% inputadata.inputdistrs=distribtouse;
% inputpars.n_distribs=numel(inputadata.inputdistrs);
% inputpars.dirstrcenters=(1:inputpars.n_distribs);
% inputpars.xtickslabelvector={'component','pattern','unclassified'};
% inputpars.distrcolors{1}=[50,200,0]./255;
% inputpars.distrcolors{2}=[255,150,0]./255;
% inputaxh = plot_violinplot_PN_new(inputadata,inputpars);
% [~,pval_csi_pc] = ttest2(distribtouse{1}',distribtouse{2}');
% [~,pval_csi_c] = ttest(distribtouse{1});
% [~,pval_csi_p] = ttest(distribtouse{2});
% [~,pval_csi_unc] = ttest(distribtouse{3});
% usedxlim=get(gca,'xlim');
% hold on;
% plot(gca,usedxlim,[0,0],'--','linewidth',2,'color',[0.5,0.5,0.5])
% text(gca,-0.25,-0.5,['comp vs. patt mean dsi diff p = ',num2str(pval_csi_pc)],'fontsize',12)
% text(gca,-0.25,-0.55,['comp mean dsi p = ',num2str(pval_csi_c)],'fontsize',12)
% text(gca,-0.25,-0.6,['patt mean dsi p = ',num2str(pval_csi_p)],'fontsize',12)
% text(gca,-0.25,-0.65,['uncl mean dsi p = ',num2str(pval_csi_unc)],'fontsize',12)
% xtickangle(45)
% set(gca,'fontsize',12)
% axis square

% plot distributions of CSI per area ---
subplot(1,2,2)
% decide wheter to use max or unrolled
distribtouse=CSI_to_use_per_area(1:3); % collapsed_max_fitted_rfs_r2_distribs
inputpars.inputaxh=gca;
hold(inputpars.inputaxh,'on')
% set settings for violin distribution plotting
inputpars.boxplotwidth=0.4;%0.5;
inputpars.boxplotlinewidth=2;
inputpars.densityplotwidth=0.5;%0.5;
inputpars.yimtouse=[-0.65,0.65];
% inputpars.yimtouse=[0,8];
inputpars.scatterjitter=inputpars.boxplotlinewidth*0.1;
inputpars.scatteralpha=0.15;
inputpars.scattersize=20;
inputpars.distralpha=0.5;
inputpars.xlabelstring=[];
inputpars.ylabelstring='CSI';
inputpars.titlestring=['dsi (',nametagbis,')',' ( comp # = ',...
    num2str(numel(distribtouse{1}),'%.0f'),...
    ' - patt #  = ',num2str(numel(distribtouse{2}),'%.0f'),...
    ' - uncl #  = ',num2str(numel(distribtouse{3}),'%.0f'),' )'];
inputpars.boolscatteron=1;
inputpars.ks_bandwidth=0.05;
inputpars.xlimtouse=[-0.5,4.5]; %[-1,5];
% plot violins
inputadata.inputdistrs=distribtouse;
inputpars.n_distribs=numel(inputadata.inputdistrs);
inputpars.dirstrcenters=(1:inputpars.n_distribs);
inputpars.xtickslabelvector={'V1','LM','RL'};
inputpars.distrcolors{1}=[255,255,255].*(0.75)./255;
inputpars.distrcolors{2}=[255,255,255].*(0.55)./255;
inputpars.distrcolors{2}=[255,255,255].*(0.35)./255;
inputaxh = plot_violinplot_PN_new(inputadata,inputpars); %#ok<NASGU>
[~,pval_csi_V1_LM] = ttest2(distribtouse{1}',distribtouse{2}');
[~,pval_csi_V1_RL] = ttest2(distribtouse{1}',distribtouse{3}');
[~,pval_csi_RL_LM] = ttest2(distribtouse{2}',distribtouse{3}');
usedxlim=get(gca,'xlim');
hold on;
plot(gca,usedxlim,[0,0],'--','linewidth',2,'color',[0.5,0.5,0.5])
text(gca,-0.25,-0.5,['V1 LM mean dsi diff p = ',num2str(pval_csi_V1_LM)],'fontsize',12)
text(gca,-0.25,-0.45,['V1 RL mean dsi diff p = ',num2str(pval_csi_V1_RL)],'fontsize',12)
text(gca,-0.25,-0.40,['RL LM mean dsi diff p = ',num2str(pval_csi_RL_LM)],'fontsize',12)
xtickangle(45)
set(gca,'fontsize',12)
axis square

saveas(fighand0,[outfold,filesep,nametag,'_scatter'],'jpg')
print(fighand0,'-depsc','-painters',[outfold,filesep,nametag,'_scatter','.eps'])
close all

%% perform rdm analysis and store representations -------------------------
% NB:revision additions from here on (26/06/2023)

addpath E:\Backups\Personal_bk\DorsalNet

% initialize RDM analyis inputs
obs_COUNT_per_class=cell(size(collapsed_obs_COUNT_distribs));
pred_COUNT_per_class=cell(size(collapsed_pred_COUNT_distribs));
% get number of stim conditions
ndir=size(collapsed_obs_COUNT_distribs{1},1);
ntype=size(collapsed_obs_COUNT_distribs{1},3);
% loop over cell classes
for cell_class_idx=1:numel(cell_types_codes)
    % get permuted current count mat - observed
    current_countmat=permute(collapsed_obs_COUNT_distribs{cell_class_idx},[2,1,3]);
    % store permuted current count mat - observed
    obs_COUNT_per_class{cell_class_idx}=normailze_along_dim(reshape(current_countmat,[size(current_countmat,1),ndir*ntype]),1);
    % get permuted current count mat - observed
    current_countmat=permute(collapsed_pred_COUNT_distribs{cell_class_idx},[2,1,3]);
    % store permuted current count mat - observed
    pred_COUNT_per_class{cell_class_idx}=normailze_along_dim(reshape(current_countmat,[size(current_countmat,1),ndir*ntype]),1);
    % store n neu per class
    n_per_class(cell_class_idx)=size(current_countmat,1); %#ok<SAGROW>
end

% set smoothpar for represenattion preprocessing smoothpar
smoothpar=0.5;
smoothbool=1;

% get matrix of median difference surprises
obs_RDM_per_class=cell(size(1,numel(n_per_class)));
pred_RDM_per_class=cell(size(1,numel(n_per_class)));
for class_idx=1:numel(n_per_class)
    % initialize rdm for current class
    obs_RDM_per_class{class_idx}=NaN(ndir*ntype,ndir*ntype);
    pred_RDM_per_class{class_idx}=NaN(ndir*ntype,ndir*ntype);
    % get preprocessed representations
    if smoothbool
        obs_rep_mat=max_normalize_halves(gaussianSmooth1D(obs_COUNT_per_class{class_idx}, smoothpar, 1));
        pred_rep_mat=max_normalize_halves(gaussianSmooth1D(pred_COUNT_per_class{class_idx}, smoothpar, 1));
    else
        obs_rep_mat=obs_COUNT_per_class{class_idx};
        pred_rep_mat=pred_COUNT_per_class{class_idx};
    end
    for cond1_idx=1:ndir*ntype
        for cond2_idx=1:ndir*ntype
            % compute current rdm element for current class - observed
            curr_obs_rep1=obs_rep_mat(:,cond1_idx);
            curr_obs_rep2=obs_rep_mat(:,cond2_idx);
            curr_corrval=corr(curr_obs_rep1,curr_obs_rep2);
            obs_RDM_per_class{class_idx}(cond1_idx,cond2_idx)=1-curr_corrval;
            % compute current rdm element for current class - predicted
            curr_pred_rep1=pred_rep_mat(:,cond1_idx);
            curr_pred_rep2=pred_rep_mat(:,cond2_idx);
            curr_corrval=corr(curr_pred_rep1,curr_pred_rep2);
            pred_RDM_per_class{class_idx}(cond1_idx,cond2_idx)=1-curr_corrval;
        end
    end
end

% prepare lablels for plotting rdms
labels = cell(ndir,ntype);
ordertags = NaN(ndir,ntype);
for i = 1:12
    for j = 1:2
        number = num2str((i-1)*30);
        if j == 1
            letter = 'g';
        else
            letter = 'p';
        end
        labels{i,j} = [number, letter];
        ordertags(i,j) = i+0.5*(j-1);
    end
end
labels=labels(:);
% get ideal (i.e. pattern) reordering permutation and block labels
ordertags=ordertags(:);
[~,ideal_sorting_perm] = sort(ordertags);
classlabels={'component','pattern','unclassified'};
blocklabels=floor(ordertags);
% prepare colors
inputpars.distrcolors{1}=[50,200,0]./255;
inputpars.distrcolors{2}=[255,150,0]./255;
inputpars.distrcolors{3}=[0,0,0]./(3*255);

% inpect pattern and component representations - observed
fighand3=figure('units','normalized','outerposition',[0 0 1 1]);
vartypes={'obs comp','obs patt'};
for class_idx=1:(numel(n_per_class)-1)
    % get current imputs
    mat_to_use=obs_RDM_per_class{class_idx};
    colortouse=inputpars.distrcolors{class_idx};
    vartype=vartypes{class_idx};
    % compute block modularity - pattern
    current_blk_mod = compute_block_modularity_PN(mat_to_use, blocklabels);
    % plot dendrogram
    subplot(2,2,2+(class_idx-1)*2)
    var_Z = linkage(mat_to_use,'average','euclidean');
    % perform hierarchical clustering
    [h,~,var_outperm] = dendrogram(var_Z);
    % compute block modularity - best
    dendrogram_height_cutoff = 0.9;
    var_T = cluster(var_Z, 'cutoff', dendrogram_height_cutoff, 'criterion', 'distance');
    current_blk_mod_best = compute_block_modularity_PN(mat_to_use, var_T);
    hold on
    plot([0,ndir*ntype],[dendrogram_height_cutoff,dendrogram_height_cutoff],'--','Linewidth',2,'Color',colortouse)
    xticklabels(gca,labels(var_outperm))
    set( h, 'Color', 'k' );
    set( h, 'LineWidth', 3 );
    set(gca,'fontsize',12)
    xlim(gca,[0,ndir*ntype]);
    if smoothbool
        ylim(gca,[0,2.5]);
    else
        ylim(gca,[0,2]);
    end
    xtickangle(gca,30)
    title([vartype,' - HC dendrogram - b mod = ',num2str(round(current_blk_mod_best,3))],'Color',colortouse)
    % plot rdm
    sorting_perm_to_use=ideal_sorting_perm;
    subplot(2,2,1+(class_idx-1)*2)
    mat_to_use_reord=mat_to_use(sorting_perm_to_use,:);
    mat_to_use_reord=mat_to_use_reord(:,sorting_perm_to_use);
    overall_labels_reord=labels(sorting_perm_to_use);
    hold on;
    im1=imagesc(flipud(mat_to_use_reord)); colormap('gray');
    axx=get(im1).Parent;
    if smoothbool
        caxis(gca,[0,1.5]);
    else
        caxis(gca,[0,1])
    end
    set(axx,'YTick',1:numel(overall_labels_reord));
    set(axx,'XTick',1:numel(overall_labels_reord));
    xlim(axx,[1-0.5,numel(overall_labels_reord)+0.5]);
    ylim(axx,[1-0.5,numel(overall_labels_reord)+0.5]);
    set(axx,'YTickLabel',flipud(overall_labels_reord));
    set(axx,'XTickLabel',overall_labels_reord);
    cb1=colorbar;
    ylabel(cb1,'dissimilarity (1-corr)')
    set(cb1,'fontsize',14)
    xtickangle(90)
    set(axx,'fontsize',12)
    axis square
    title([vartype,' - HC-reordered RDM - p mod = ',num2str(round(current_blk_mod,3))],'Color',colortouse)
end
suptitle('Rat - patt and comp observed representations - all layers')
saveas(fighand3,[outfold,filesep,'Rat_obs_rdm'],'jpg')
print(fighand3,'-depsc','-painters',[[outfold,filesep,'Rat_obs_rdm'],'.eps'])
% close all

% inpect pattern and component representations - predicted
fighand4=figure('units','normalized','outerposition',[0 0 1 1]);
vartypes={'pred comp','pred patt'};
% blk_mod_pred=cell(1,(numel(n_per_class)-1));
% blklbl_pred=cell(1,(numel(n_per_class)-1));
% blk_mod_best_pred=cell(1,(numel(n_per_class)-1));
% blklbl_best_pred=cell(1,(numel(n_per_class)-1));
for class_idx=1:(numel(n_per_class)-1)
    % get current imputs
    mat_to_use=pred_RDM_per_class{class_idx};
    colortouse=inputpars.distrcolors{class_idx};
    vartype=vartypes{class_idx};
    % compute block modularity
    current_blk_mod = compute_block_modularity_PN(mat_to_use, blocklabels);
    %     blklbl_pred{class_idx}=blocklabels;
    %     blk_mod_pred{class_idx}=current_blk_mod;
    % plot dendrogram
    subplot(2,2,2+(class_idx-1)*2)
    var_Z = linkage(mat_to_use,'average','euclidean');
    % perform hierarchical clustering
    [h,~,var_outperm] = dendrogram(var_Z);
    % compute block modularity - best
    dendrogram_height_cutoff = 1;
    var_T = cluster(var_Z, 'cutoff', dendrogram_height_cutoff, 'criterion', 'distance');
    current_blk_mod_best = compute_block_modularity_PN(mat_to_use, var_T);
    %     blklbl_best_pred{class_idx}=var_T;
    %     blk_mod_best_pred{class_idx}=current_blk_mod_best;
    hold on
    plot([0,ndir*ntype],[dendrogram_height_cutoff,dendrogram_height_cutoff],'--','Linewidth',2,'Color',colortouse)
    xticklabels(gca,labels(var_outperm))
    set( h, 'Color', 'k' );
    set( h, 'LineWidth', 3 );
    set(gca,'fontsize',12)
    xlim(gca,[0,ndir*ntype]);
    if smoothbool
        ylim(gca,[0,2.5]);
    else
        ylim(gca,[0,2]);
    end
    xtickangle(gca,30)
    title([vartype,' - HC dendrogram - b mod = ',num2str(round(current_blk_mod_best,3))],'Color',colortouse)
    % plot rdm
    sorting_perm_to_use=ideal_sorting_perm;
    subplot(2,2,1+(class_idx-1)*2)
    mat_to_use_reord=mat_to_use(sorting_perm_to_use,:);
    mat_to_use_reord=mat_to_use_reord(:,sorting_perm_to_use);
    overall_labels_reord=labels(sorting_perm_to_use);
    hold on;
    im1=imagesc(flipud(mat_to_use_reord)); colormap('gray');
    axx=get(im1).Parent;
    if smoothbool
        caxis(gca,[0,1.5]);
    else
        caxis(gca,[0,1])
    end
    set(axx,'YTick',1:numel(overall_labels_reord));
    set(axx,'XTick',1:numel(overall_labels_reord));
    xlim(axx,[1-0.5,numel(overall_labels_reord)+0.5]);
    ylim(axx,[1-0.5,numel(overall_labels_reord)+0.5]);
    set(axx,'YTickLabel',flipud(overall_labels_reord));
    set(axx,'XTickLabel',overall_labels_reord);
    cb1=colorbar;
    ylabel(cb1,'dissimilarity (1-corr)')
    set(cb1,'fontsize',14)
    xtickangle(90)
    set(axx,'fontsize',12)
    axis square
    title([vartype,' - HC-reordered RDM - p mod = ',num2str(round(current_blk_mod,3))],'Color',colortouse)
end
suptitle('Rat patt and comp LN-predicted representations - all layers')
saveas(fighand4,[outfold,filesep,'Rat_pred_rdm'],'jpg')
print(fighand4,'-depsc','-painters',[[outfold,filesep,'Rat_pred_rdm'],'.eps'])
close all

%store bullet distr data
bullet_distr.idx_consist=idx_consist;
bullet_distr.collapsed_Zc_O_distribs=collapsed_Zc_O_distribs;
bullet_distr.collapsed_Zp_O_distribs=collapsed_Zp_O_distribs;
bullet_distr.collapsed_Zc_P_distribs=collapsed_Zc_P_distribs;
bullet_distr.collapsed_Zp_P_distribs=collapsed_Zp_P_distribs;

% save results of rdm analysis
save([outfold,filesep,'RDM_datastructure_Rat.mat'],...
    'obs_COUNT_per_class',...
    'pred_COUNT_per_class',...
    'obs_RDM_per_class',...
    'pred_RDM_per_class',...
    'labels',...
    'classlabels',...
    'blocklabels',...
    'ideal_sorting_perm',...
    'bullet_distr');

% NB (26/06/2023): this datastructure is to be used
% "perform_representation_regression_DorsalNet_PN.m" in the other branch
% after displacing the .mat file in the other folder
% "E:\Backups\Personal_bk\DorsalNet"

%% analyze plaid & grating modulation indeces and produce raster examples ----------------------

% add path to recent (Mattia's) modulation index code
addpath('D:\Backups\Personal_bk\PN_acute_analysis\scripts\new_analyses_2022\revision_MI_code')
cell_class_labels={'component','pattern','unclassified'};
% set whether to rerun computation of MI
bool_rerun_mi_analysis=0;
% set whether to plot single neuron diagnostics
bool_plot_mi_diagnostics=0;

if not(bool_rerun_mi_analysis) && exist([outfold,filesep,'MI_results_',labelstouse,'.mat']) %#ok<EXIST>
    
    % load previously computaed MI
    load([outfold,filesep,'MI_results_',labelstouse,'.mat'])
    
elseif bool_rerun_mi_analysis
    
    % load original tuning data
    tuning_file=['D:\Backups\Personal_bk\PN_acute_analysis\processed_data',filesep,'Tuning','.mat']; %#ok<*UNRCH>
    tuning=load(tuning_file);
    
    % get parameters
    sr=videoframerate;
    countingwindowlim=[0,1];
    input_dirs=DIR;
    
    % set interp factor
    input_interp_factor=interp_factor;
    
    % set best spike count and CSI structures
    mod_MIf1f0=cell(1,numel(cell_types_codes));
    mod_MIf1z=cell(1,numel(cell_types_codes));
    mod_spectrum=cell(1,numel(cell_types_codes));
    mod_FR=cell(1,numel(cell_types_codes));
    mod_COUNT=cell(1,numel(cell_types_codes));
    mod_RASTER=cell(1,numel(cell_types_codes));
    mod_pDIR=cell(1,numel(cell_types_codes));
    mod_SF=cell(1,numel(cell_types_codes));
    mod_TF=cell(1,numel(cell_types_codes));
    mod_nNUM=cell(1,numel(cell_types_codes));
    
    % loop over cell classes
    for cell_class_idx=1:numel(cell_types_codes)
        
        % get neurons id
        current_ns=collapsed_ns_distribs{cell_class_idx};
        % get SF and TF (same of grating best condition)
        current_g_pSFs=pref_SF(current_ns,1);
        current_g_pTFs=pref_TF(current_ns,1);
        current_p_pSFs=pref_SF(current_ns,1);
        current_p_pTFs=pref_TF(current_ns,1);
        
        % initialize output matrices
        mod_MIf1f0{cell_class_idx}=cell(1,numel(current_ns));
        mod_MIf1z{cell_class_idx}=cell(1,numel(current_ns));
        mod_spectrum{cell_class_idx}=cell(1,numel(current_ns));
        mod_FR{cell_class_idx}=cell(1,numel(current_ns));
        mod_COUNT{cell_class_idx}=cell(1,numel(current_ns));
        mod_RASTER{cell_class_idx}=cell(1,numel(current_ns));
        mod_pDIR{cell_class_idx}=cell(1,numel(current_ns));
        mod_SF{cell_class_idx}=cell(1,numel(current_ns));
        mod_TF{cell_class_idx}=cell(1,numel(current_ns));
        mod_nNUM{cell_class_idx}=cell(1,numel(current_ns));
        
        % loop over neurons
        for i=1:numel(current_ns)
            
            % initialize storage structures
            mod_MIf1f0{cell_class_idx}{i}=NaN(length(input_dirs),numel(stimulustypes));
            mod_MIf1z{cell_class_idx}{i}=NaN(length(input_dirs),numel(stimulustypes));
            mod_spectrum{cell_class_idx}{i}=cell(length(input_dirs),numel(stimulustypes));
            mod_FR{cell_class_idx}{i}=NaN(sr*diff(countingwindowlim),length(input_dirs),numel(stimulustypes));
            mod_COUNT{cell_class_idx}{i}=NaN(1,length(input_dirs),numel(stimulustypes));
            mod_RASTER{cell_class_idx}{i}=cell(1,length(input_dirs),numel(stimulustypes));
            mod_pDIR{cell_class_idx}{i}=cell(1,numel(stimulustypes));
            
            tic
            % loop over stim types
            for k=1:numel(stimulustypes)
                % get current cell class
                current_label=cell_class_idx;
                % get current sf and tf to use
                if k==1
                    input_sf=current_g_pSFs(i);
                    input_tf=current_g_pTFs(i);
                elseif k==2
                    input_sf=current_p_pSFs(i);
                    input_tf=current_p_pTFs(i);
                end
                input_neuron_num=current_ns(i);
                input_dirs=DIR;
                input_stimtype=stimulustypes{k};
                % loop over directions
                for dir_idx=1:length(input_dirs)
                    % get current direction
                    current_direction=input_dirs(dir_idx);
                    % get observed psth
                    [current_psth,input_RASTER,current_psth_edges] = get_psth_PN(input_neuron_num,input_sf,input_tf,current_direction,input_stimtype);
                    % get counting window span
                    counting_window_span=countingwindowlim(2)-countingwindowlim(1);
                    % smooth observed psth to match timescale of prediction
                    gaussFilter = gausswin(10,1/0.5);
                    gaussFilter = gaussFilter / sum(gaussFilter);
                    current_psth=conv(current_psth, gaussFilter,'same');
                    % resample observed psth
                    current_psth=squeeze(current_psth);
                    new_samplenum=sr.*(counting_window_span);
                    old_samplenum=length(current_psth);
                    current_psth_resampled = resample(current_psth,new_samplenum,old_samplenum);
                    % compute phase modulation
                    psth=current_psth_resampled;
                    psth_time=(0:length(current_psth_resampled)).*(1/sr);
                    tf_target=pref_TF(current_ns(i),k);
                    boolplot=0;
                    [ current_MIf1z, current_MIf1f0 , current_MIplothandle, current_MIplotdata ] = ...
                        get_modulation_index( psth, max(psth_time), tf_target, boolplot );
                    % store results
                    mod_FR{cell_class_idx}{i}(:,dir_idx,k)=current_psth_resampled;
                    mod_COUNT{cell_class_idx}{i}(dir_idx,k)=nansum(current_psth_resampled);
                    mod_RASTER{cell_class_idx}{i}{dir_idx,k}=input_RASTER;
                    mod_MIf1f0{cell_class_idx}{i}(dir_idx,k)=current_MIf1f0;
                    mod_MIf1z{cell_class_idx}{i}(dir_idx,k)=current_MIf1z;
                    mod_spectrum{cell_class_idx}{i}{dir_idx,k}=current_MIplotdata;
                    mod_SF{cell_class_idx}{i}=input_sf;
                    mod_TF{cell_class_idx}{i}=input_tf;
                    mod_nNUM{cell_class_idx}{i}=input_neuron_num;
                end
                % compute and store direction
                mod_pDIR{cell_class_idx}{i}{k}=DIR(find(mod_COUNT{cell_class_idx}{i}(:,k)==max(mod_COUNT{cell_class_idx}{i}(:,k))));
            end
            toc
            
            if bool_plot_mi_diagnostics
                % initialize figure
                ffffffff = figure('units','normalized','outerposition',[0 0 1 1]);
                % plot psth and raster ----------------------------
                axis square
                % get data to plot for current neuron
                curr_pDIR=mod_pDIR{cell_class_idx}{i};
                curr_SF=mod_SF{cell_class_idx}{i};
                curr_TF=mod_TF{cell_class_idx}{i};
                curr_nNUM=mod_nNUM{cell_class_idx}{i};
                curr_FR=mod_FR{cell_class_idx}{i};
                curr_RASTER=mod_RASTER{cell_class_idx}{i};
                current_MIplotdata=mod_spectrum{cell_class_idx}{i};
                curr_COUNT=mod_COUNT{cell_class_idx}{i};
                % draw raster
                for k=1:2
                    % draw psth ------------
                    sb1=subplot(2,3,2+(k-1));
                    hold on;
                    psth_observed=curr_FR(:,DIR==curr_pDIR{k},k);
                    psth_observed=psth_observed./(max(psth_observed(:)));
                    % set color and tag to use
                    if current_label==2
                        coltuse=[255,150,0]./255;
                    elseif current_label==1
                        coltuse=[50,200,0]./255;
                    elseif current_label==3
                        coltuse=[150,150,150]./255;
                    end
                    hold on;
                    % draw psth
                    pipi=plot(gca,(0:length(psth_observed)).*(1/sr),[psth_observed(1);psth_observed],'-','Color',coltuse./(k),'LineWidth',2.5);
                    plot_shaded_auc(gca,(0:length(psth_observed))*(1/sr),[psth_observed(1);psth_observed]',0.15,coltuse./(k))
                    % draw raster ------------
                    raster=curr_RASTER{curr_pDIR{k}==DIR,k};
                    spp=raster;
                    spk=0;
                    for kkk=1:numel(raster)
                        if not(isempty(raster{kkk}))
                            plot(raster{kkk},1.05+0.05*kkk,'.','color',coltuse./(k),'MarkerSize',15)
                            spk=spk+numel(raster{kkk});
                        else
                        end
                    end
                    xlim([-0.2,1.2]);
                    ylim([-0,4]);
                    if k==1
                        tt=text(0.05,3.4,['pDIR grating = ',num2str(curr_pDIR{1}),' d'],'FontSize',12);
                        ttt=text(0.05,3.2,['spike count grating = ',num2str(spk)],'FontSize',12);
                    else
                        tt2=text(0.05,3.4,['pDIR plaid = ',num2str(curr_pDIR{2}),' d'],'FontSize',12);
                        ttt2=text(0.05,3.2,['spike count plaid = ',num2str(spk)],'FontSize',12);
                    end
                    plot([0,0],[0,5],'--k', 'LineWidth',2)
                    plot([1,1],[0,5],'--k', 'LineWidth',2)
                    hlabelx=get(gca,'Xlabel');
                    set(hlabelx,'String','time (s)','FontSize',12,'color','k')
                    hlabely=get(gca,'Ylabel');
                    set(hlabely,'String','normalized firing rate','FontSize',12,'color','k')
                    legend(gca,pipi,{'grating','plaid'})
                    title(['raster and psth (n=',num2str(curr_nNUM),') - DIR=','best',' - TF=',num2str(curr_TF),' - SF=',num2str(curr_SF)])
                    set(gca,'FontSize',12);
                end
                % plot polar plot tuning curve ----------------------------
                for k=1:2
                    subplot(2,3,5+(k-1));
                    % fetch data
                    N_f=current_MIplotdata{curr_pDIR{k}==DIR,k}.N_f;
                    f=current_MIplotdata{curr_pDIR{k}==DIR,k}.f;
                    fidx=current_MIplotdata{curr_pDIR{k}==DIR,k}.fidx;
                    pow=current_MIplotdata{curr_pDIR{k}==DIR,k}.pow;
                    meanspect=current_MIplotdata{curr_pDIR{k}==DIR,k}.meanspect;
                    sigspect=current_MIplotdata{curr_pDIR{k}==DIR,k}.sigspect;
                    tTF=current_MIplotdata{curr_pDIR{k}==DIR,k}.TF;
                    F1F0=current_MIplotdata{curr_pDIR{k}==DIR,k}.F1F0;
                    F1z=current_MIplotdata{curr_pDIR{k}==DIR,k}.F1z;
                    % draw spectrum
                    plot(f(1:round(N_f/2)+1),pow(1:round(N_f/2)+1),'color',coltuse./(k),'LineWidth',3);
                    hold on
                    plot(f(fidx),pow(fidx),'o','LineWidth',6,'color',0*coltuse./(k));
                    plot(f(1:round(N_f/2)+1),(meanspect+sigspect)*ones(size(pow(1:round(N_f/2)+1))),'--','color',[0.5,0.5,0.5],'LineWidth',1);
                    plot(f(1:round(N_f/2)+1),(meanspect-sigspect)*ones(size(pow(1:round(N_f/2)+1))),'--','color',[0.5,0.5,0.5],'LineWidth',1);
                    plot(f(1:round(N_f/2)+1),(meanspect)*ones(size(pow(1:round(N_f/2)+1))),'.-','color',[0.5,0.5,0.5],'LineWidth',0.5);
                    ylimit=get(gca,'ylim');
                    xlimit=get(gca,'xlim');
                    ttt=text(0.5*xlimit(2),0.85*ylimit(2),['target TF = ',num2str(tTF),' Hz'],'FontSize',12);
                    set(gca,'FontSize',10);
                    hlabelx=get(gca,'Xlabel');
                    set(hlabelx,'String','f [Hz]','FontSize',10,'color','k')
                    hlabely=get(gca,'Ylabel');
                    set(hlabely,'String','PSD','FontSize',10,'color','k')
                    title(['Power spectrum (F1z = ',num2str(F1z),', F1F0 = ',num2str(F1F0),')']);
                    hold off
                    axis square
                    set(gca,'fontsize',12);
                end
                % plot polar plot tuning curve ----------------------------
                ppol=polaraxes('Position',[-0.05,0.25,.5,.5]);
                hold on;
                for k=1:2
                    temp_obs_tc=squeeze(curr_COUNT(:,k))'./max(squeeze(curr_COUNT(:,k)));
                    obs_tc=[temp_obs_tc,temp_obs_tc(1)];
                    % draw plar plots
                    p2=polarplot(ppol,[deg2rad(DIR),2*pi],obs_tc,'-');
                    set(p2,'color',coltuse./k)
                    set(p2, 'linewidth', 3.5);
                    % draw polarscatter of pref dir
                    ps2=polarscatter(ppol,deg2rad(curr_pDIR{k}),obs_tc(find(DIR==curr_pDIR{k})),...
                        125,'markerfacecolor',coltuse./k,'markeredgecolor',coltuse./k);
                end
                title(ppol,['tuning polar plots (n=',num2str(curr_nNUM),') - DIR=',num2str(curr_pDIR{1}),' - TF=',num2str(curr_TF),' - SF=',num2str(curr_SF)])
                set(ppol,'fontsize',12);
                sgtitle(['MI diagnostics n',num2str(curr_nNUM),' - ',num2str(cell_class_labels{cell_class_idx})])
                % save figure
                saveas(ffffffff,[outfold,filesep,'MI_diagnostics_nNUM',num2str(curr_nNUM),'_n',num2str(i)],'jpg')
                print(ffffffff,'-depsc','-painters',[[outfold,filesep,'MI_diagnostics_n',num2str(curr_nNUM)],'.eps'])
                close all
            end
            
            % output advancement message
            fprintf(['neuron#',num2str(i),' class#',num2str(cell_class_idx),' analyzed...\n'])
        end
    end
    save([outfold,filesep,'MI_results_',labelstouse,'.mat'],...
        'mod_MIf1f0',...
        'mod_MIf1z',...
        'mod_spectrum',...
        'mod_FR',...
        'mod_COUNT',...
        'mod_RASTER',...
        'mod_pDIR',...
        'mod_SF',...
        'mod_TF',...
        'mod_nNUM')
end

% initalize output datastructure
F1z_distr_grating=cell(1,numel(cell_types_codes));
F1z_distr_plaid=cell(1,numel(cell_types_codes));
PI_distr=cell(1,numel(cell_types_codes));

% loop over cell classes
for cell_class_idx=1:numel(cell_types_codes)
    F1z_distr_grating{cell_class_idx}=NaN(1,numel(mod_MIf1z{cell_class_idx}));
    F1z_distr_plaid{cell_class_idx}=NaN(1,numel(mod_MIf1z{cell_class_idx}));
    PI_distr{cell_class_idx}=NaN(1,numel(mod_MIf1z{cell_class_idx}));
    % loop over neurons
    for i=1:numel(mod_MIf1z{cell_class_idx})
        F1z_distr_grating{cell_class_idx}(i)=mod_MIf1z{cell_class_idx}{i}(mod_pDIR{cell_class_idx}{i}{1}==DIR,1);
        F1z_distr_plaid{cell_class_idx}(i)=mod_MIf1z{cell_class_idx}{i}(mod_pDIR{cell_class_idx}{i}{2}==DIR,2);
    end
    % get pattern index
    PI_distr{cell_class_idx}=collapsed_Zp_O_oldbest_distribs{cell_class_idx}-collapsed_Zc_O_oldbest_distribs{cell_class_idx};
end

% loop over cell classes
simple_frac=NaN(1,numel(cell_types_codes));
simple_num=NaN(1,numel(cell_types_codes));
complex_num=NaN(1,numel(cell_types_codes));
for cell_class_idx=1:numel(cell_types_codes)
    % count simple and complex cells
    n_simple=sum(F1z_distr_grating{cell_class_idx}>3);
    n_complex=sum(F1z_distr_grating{cell_class_idx}<=3);
    n_tot=n_simple+n_complex;
    simple_frac(cell_class_idx)=n_simple./n_tot;
    simple_num(cell_class_idx)=n_simple;
    complex_num(cell_class_idx)=n_complex;
end
% run Fisher exact test
x = table(simple_num(1:2)',complex_num(1:2)','VariableNames',{'Simple','Complex'},'RowNames',{'Component','Pattern'});
[h,pfishexact] = fishertest(x);

% plot MI analysis result (violins) ---------------------------------------
fighand000=figure('units','normalized','outerposition',[0 0 1 1]);
% plot distributions of MI --- (gratings)
subplot(1,3,1)
% decide wheter to use max or unrolled
distribtouse=F1z_distr_grating(1:2);
inputpars.inputaxh=gca;
hold(inputpars.inputaxh,'on')
% set settings for violin distribution plotting
inputpars.boxplotwidth=0.4;%0.5;
inputpars.boxplotlinewidth=2;
inputpars.densityplotwidth=0.5;%0.5;
inputpars.yimtouse=[-0.5,5];
% inputpars.yimtouse=[0,8];
inputpars.scatterjitter=inputpars.boxplotlinewidth*0.1;
inputpars.scatteralpha=0.15;
inputpars.scattersize=20;
inputpars.distralpha=0.5;
inputpars.xlabelstring=[];
inputpars.ylabelstring='MI';
inputpars.titlestring=['modulation index (','MI',') - ','gratings',' ( comp # = ',...
    num2str(numel(distribtouse{1}),'%.0f'),...
    ' - patt #  = ',num2str(numel(distribtouse{2}),'%.0f'),')'];
inputpars.boolscatteron=1;
inputpars.ks_bandwidth=0.175;
inputpars.xlimtouse=[-0.5,3.5]; %[-1,5];
% plot violins
inputadata.inputdistrs=distribtouse;
inputpars.n_distribs=numel(inputadata.inputdistrs);
inputpars.dirstrcenters=(1:inputpars.n_distribs);
inputpars.xtickslabelvector={'component','pattern'};
inputpars.distrcolors{1}=[50,200,0]./255;
inputpars.distrcolors{2}=[255,150,0]./255;
inputaxh = plot_violinplot_PN_new(inputadata,inputpars); %#ok<NASGU>
hold on;
plot(get(gca,'xlim'),[3,3],'--','linewidth',2,'color',[0.5,0.5,0.5])
[pval_csi_pc,~] = ranksum(distribtouse{1}',distribtouse{2}');
median_MI_comp=nanmedian(distribtouse{1});
median_MI_patt=nanmedian(distribtouse{2});
usedxlim=get(gca,'xlim'); %#ok<NASGU>
hold on;
text(gca,0,4.25,['comp vs. patt median mi diff p = ',num2str(pval_csi_pc)],'fontsize',12)
text(gca,0,4.5,['comp median val = ',num2str(median_MI_comp)],'fontsize',12)
text(gca,0,4.75,['patt median val = ',num2str(median_MI_patt)],'fontsize',12)
xtickangle(45)
set(gca,'fontsize',12)
axis square
% plot distributions of MI --- (plaids)
subplot(1,3,2)
% decide wheter to use max or unrolled
distribtouse=F1z_distr_plaid(1:2);
inputpars.inputaxh=gca;
hold(inputpars.inputaxh,'on')
% set settings for violin distribution plotting
inputpars.boxplotwidth=0.4;%0.5;
inputpars.boxplotlinewidth=2;
inputpars.densityplotwidth=0.5;%0.5;
inputpars.yimtouse=[-0.5,5];
% inputpars.yimtouse=[0,8];
inputpars.scatterjitter=inputpars.boxplotlinewidth*0.1;
inputpars.scatteralpha=0.15;
inputpars.scattersize=20;
inputpars.distralpha=0.5;
inputpars.xlabelstring=[];
inputpars.ylabelstring='MI';
inputpars.titlestring=['modulation index (','MI',') - ','plaids',' ( comp # = ',...
    num2str(numel(distribtouse{1}),'%.0f'),...
    ' - patt #  = ',num2str(numel(distribtouse{2}),'%.0f'),')'];
inputpars.boolscatteron=1;
inputpars.ks_bandwidth=0.175;
inputpars.xlimtouse=[-0.5,3.5]; %[-1,5];
% plot violins
inputadata.inputdistrs=distribtouse;
inputpars.n_distribs=numel(inputadata.inputdistrs);
inputpars.dirstrcenters=(1:inputpars.n_distribs);
inputpars.xtickslabelvector={'component','pattern'};
inputpars.distrcolors{1}=([50,200,0]./255)/2;
inputpars.distrcolors{2}=([255,150,0]./255)/2;
inputaxh = plot_violinplot_PN_new(inputadata,inputpars); %#ok<NASGU>
hold on;
plot(get(gca,'xlim'),[3,3],'--','linewidth',2,'color',[0.5,0.5,0.5])
[pval_csi_pc,~] = ranksum(distribtouse{1}',distribtouse{2}');
median_MI_comp=nanmedian(distribtouse{1});
median_MI_patt=nanmedian(distribtouse{2});
usedxlim=get(gca,'xlim');
hold on;
text(gca,0,4.25,['comp vs. patt median mi diff p = ',num2str(pval_csi_pc)],'fontsize',12)
text(gca,0,4.5,['comp median val = ',num2str(median_MI_comp)],'fontsize',12)
text(gca,0,4.75,['patt median val = ',num2str(median_MI_patt)],'fontsize',12)
xtickangle(45)
set(gca,'fontsize',12)
axis square
subplot(1,3,3)
inputpars.distrcolors{1}=[50,200,0]./255;
inputpars.distrcolors{2}=[255,150,0]./255;
totnum=simple_num+complex_num;
simplefrac=simple_num./totnum;
barposition=[1,2];
hold on;
bar(barposition(1),simplefrac(1),'facecolor',inputpars.distrcolors{1},...
    'edgecolor',inputpars.distrcolors{1})
bar(barposition(2),simplefrac(2),'facecolor',inputpars.distrcolors{2},...
    'edgecolor',inputpars.distrcolors{2})
ylabel('fraction')
xticks(barposition)
xticklabels({'comp-simp', 'patt-simp'})
hold on;
xtickangle(45)
ylabel('fraction')
set(gca,'fontsize',12)
xlim([0,3])
axis square
ylimused=get(gca,'ylim');
text(gca,1.6,0.9.*ylimused(end),['# comp-simp = ',num2str(simple_num(1))],'fontsize',12)
text(gca,1.6,0.86.*ylimused(end),['# patt-simp = ',num2str(simple_num(2))],'fontsize',12)
text(gca,1.6,0.82.*ylimused(end),['p fisher exact = ',num2str(round(pfishexact,4))],'fontsize',12)
title(['fraction of simple cells - ','gratings',' ( comp # = ',...
    num2str(numel(distribtouse{1}),'%.0f'),...
    ' - patt #  = ',num2str(numel(distribtouse{2}),'%.0f'),')']);
% add suptitile
sgtitle('grating and plaid MI distributions (pattern vs. components)')
saveas(fighand000,[outfold,filesep,'MI_distributions'],'jpg')
print(fighand000,'-depsc','-painters',[outfold,filesep,'MI_distributions','.eps'])
close all

% plot MI analysis result (PI scatter) ------------------------------------
fighand0000=figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
inputpars.distrcolors{1}=[50,200,0]./255;
inputpars.distrcolors{2}=[255,150,0]./255;
inputpars.distrcolors{3}=0.5*[255,255,255]./255;
scatter(PI_distr{3},F1z_distr_grating{3},100,'MarkerFaceColor',inputpars.distrcolors{3},'MarkerEdgeColor',inputpars.distrcolors{3})
scatter(PI_distr{1},F1z_distr_grating{1},100,'MarkerFaceColor',inputpars.distrcolors{1},'MarkerEdgeColor',inputpars.distrcolors{1})
scatter(PI_distr{2},F1z_distr_grating{2},100,'MarkerFaceColor',inputpars.distrcolors{2},'MarkerEdgeColor',inputpars.distrcolors{2})
plot([0,0],get(gca,'ylim'),'--','linewidth',2,'color',[0.5,0.5,0.5])
ylabel('MI')
xlabel('PI')
PI_distr_all=[PI_distr{1}',PI_distr{2}',PI_distr{3}']';
F1z_distr_grating_all=[F1z_distr_grating{1},F1z_distr_grating{2},F1z_distr_grating{3}]';
[coor_v,corr_p]=corr(PI_distr_all,F1z_distr_grating_all);
text(gca,2,0.4,['PI vs. MI corr v = ',num2str(round(coor_v,2))],'fontsize',12)
text(gca,2,0.2,['PI vs. MI corr p = ',num2str(corr_p)],'fontsize',12)
title(' modulation index vs. pattern index scatter (pattern vs. components) - gratings')
set(gca,'fontsize',12)
axis square
saveas(fighand0000,[outfold,filesep,'MI_PI_scatter'],'jpg')
print(fighand0000,'-depsc','-painters',[outfold,filesep,'MI_PI_scatter','.eps'])
close all

%% analyze STA best frames population average and equiv. SF distributions ----------------------

% initialize population STA structures per area
selectedframes=cell(1,numel(target_areas));
STA_bestframe=cell(1,numel(target_areas));
selectedframes_contrast=cell(1,numel(target_areas));
pop_STAs=cell(1,numel(target_areas));
pop_wSTA=cell(1,numel(target_areas));
pop_wSTA_weights=cell(1,numel(target_areas));
pop_rSTA=cell(1,numel(target_areas));
pop_rSTA_std=cell(1,numel(target_areas));
pop_classid=cell(1,numel(target_areas));
pop_TCs=cell(1,numel(target_areas));
pop_ns=cell(1,numel(target_areas));
% loop over areas
for target_area_idx=1:numel(target_areas)
    tic
    % get current sta data filepath
    current_filepath=[data_folder,filesep,'RFs_new_datastructure_',target_areas{target_area_idx}];
    % load STA data
    current_STA_data=load(current_filepath);
    % get frame selection
    current_selectedframes=current_STA_data.bestfr_contrast;
    % initialize best STA frames storage
    current_STA_bestframe=NaN(size(current_STA_data.Zwsta,1),size(current_STA_data.Zwsta,2),numel(current_selectedframes));
    % initialize contrast of best STA frames storage
    current_selectedframes_contrast=NaN(1,numel(current_selectedframes));
    % loop over neurons
    for i=1:numel(current_selectedframes)
        % get best STA frames
        current_STA_bestframe(:,:,i)=current_STA_data.Zwsta(:,:,current_selectedframes(i),i);
        % get contrast of best STA frames
        current_selectedframes_contrast(i)=current_STA_data.contrast(current_selectedframes(i),i);
    end
    % compute contrast weighted population STA
    weightingvid=permute(repmat(current_selectedframes_contrast,[size(current_STA_bestframe,1),1,size(current_STA_bestframe,2)]),[1,3,2]);
    weightingvid=weightingvid./sum(current_selectedframes_contrast);
    curr_pop_wSTA=nansum(abs(current_STA_bestframe).*weightingvid,3);
    % compute contrast raw population STA
    curr_pop_rSTA=nanmean(abs(current_STA_bestframe),3);
    % compute contrast raw population STA std
    curr_pop_rSTA_std=nanstd(abs(current_STA_bestframe),[],3);
    % get current neuron numbers
    current_ns = current_STA_data.neuron_number;
    pop_ns{target_area_idx} = current_ns;
    % initielize neuron class labels
    pop_classid{target_area_idx} = NaN(1,numel(current_selectedframes));
    % fetch STAs amd tuning curves for current population
    pop_TCs{target_area_idx} = tuning_curve(:,current_ns,:,:,:);
    pop_STAs{target_area_idx} = current_STA_bestframe;
    % initialize neuron class labels
    for cell_type_idx_inner = 1:numel(cell_types_codes)
        current_cell_type_ns = ns_distribs{cell_type_idx_inner,target_area_idx};
        for ii=1:numel(current_cell_type_ns)
            %store albel for current neuron
            current_idx=find(current_cell_type_ns(ii)==current_ns);
            pop_classid{target_area_idx}(current_idx)=cell_type_idx_inner;
        end
    end
    % store results for current area
    selectedframes{target_area_idx}=current_selectedframes;
    STA_bestframe{target_area_idx}=current_STA_bestframe;
    selectedframes_contrast{target_area_idx}=current_selectedframes_contrast;
    pop_wSTA{target_area_idx}=curr_pop_wSTA;
    pop_wSTA_weights{target_area_idx}=weightingvid;
    pop_rSTA{target_area_idx}=curr_pop_rSTA;
    pop_rSTA_std{target_area_idx}=curr_pop_rSTA_std;
    toc
end
% initialize fit parameters structure
fitParams_pop_wSTA=cell(1,numel(target_areas));
fitParams_pop_rSTA=cell(1,numel(target_areas));
fittedGaussian_pop_wSTA=cell(1,numel(target_areas));
fittedGaussian_pop_rSTA=cell(1,numel(target_areas));
% loop over areas
for target_area_idx=1:numel(target_areas)
    regpars=1*ones(1,5);
    inputdata=pop_wSTA{target_area_idx};
    [ fitParams_pop_wSTA{target_area_idx}, fittedGaussian_pop_wSTA{target_area_idx}, ~ ] =...
        get_2D_gaussian_fit(inputdata-nanmean(inputdata(:)),regpars);
    inputdata=pop_rSTA{target_area_idx};
    [ fitParams_pop_rSTA{target_area_idx}, fittedGaussian_pop_rSTA{target_area_idx}, bfCost ] =...
        get_2D_gaussian_fit(inputdata-nanmean(inputdata(:)),regpars);
end
% initialize single neuron fitting results structures
goodcontrast_neuid=cell(1,numel(target_areas));
fitParams_bestSTAlobe=cell(1,numel(target_areas));
fittedGaussian_bestSTAlobe=cell(1,numel(target_areas));
SF_equiv_bestSTAlobe=cell(1,numel(target_areas));
X_center_bestSTAlobe=cell(1,numel(target_areas));
Y_center_bestSTAlobe=cell(1,numel(target_areas));
classid_bestSTAlobe=cell(1,numel(target_areas));
% set contrast threshold to use for fitting
contrast_th_to_use=7.5;
% set whether to plot fitting diagnostics
boolplotfitdiagnostics=0;
% set regularization to use for fitting
regpars=0*ones(1,5);
% loop over areas
for target_area_idx=1:numel(target_areas)
    tic
    % select good contrast STA neuron indeces
    goodcontrast_neuid{target_area_idx}=find(selectedframes_contrast{target_area_idx}>=contrast_th_to_use);
    % initialize current fit datastructures
    fitParams_bestSTAlobe{target_area_idx}=cell(1,numel(goodcontrast_neuid{target_area_idx}));
    fittedGaussian_bestSTAlobe{target_area_idx}=cell(1,numel(goodcontrast_neuid{target_area_idx}));
    SF_equiv_bestSTAlobe{target_area_idx}=NaN(1,numel(goodcontrast_neuid{target_area_idx}));
    X_center_bestSTAlobe{target_area_idx}=NaN(1,numel(goodcontrast_neuid{target_area_idx}));
    Y_center_bestSTAlobe{target_area_idx}=NaN(1,numel(goodcontrast_neuid{target_area_idx}));
    classid_bestSTAlobe{target_area_idx}=NaN(1,numel(goodcontrast_neuid{target_area_idx}));
    % loop over nezrons with good STA in current area
    for neu_idx=1:numel(goodcontrast_neuid{target_area_idx})
        % get current neuron
        current_neuid=goodcontrast_neuid{target_area_idx}(neu_idx);
        % gest best STA frame for current neuron
        inputdata_temp=imresize(STA_bestframe{target_area_idx}(:,:,current_neuid),2);
        inputdata=inputdata_temp;
        % isolate strongest lobe
        [tempmaxval,tempid]=max([max(inputdata(:)),abs(min(inputdata(:)))]);
        if tempid==1 % excitatory lobe is stronger
            inputdata(inputdata<=0)=0;
        elseif tempid==2 % inhibitory lobe is stronger
            inputdata(inputdata>=0)=0;
        end
        inputdata=abs(inputdata);
        % perform fitting
        [fitParams_bestSTAlobe{target_area_idx}{neu_idx}, fittedGaussian_bestSTAlobe{target_area_idx}{neu_idx}, ~] = get_2D_gaussian_fit(inputdata-nanmean(inputdata(:)),regpars);
        % extract Rf measure from fit
        measure1=fitParams_bestSTAlobe{target_area_idx}{neu_idx}(3);
        measure2=fitParams_bestSTAlobe{target_area_idx}{neu_idx}(4);
        % compute average bump sigma from fit
        sigmaequiv=(measure1+measure2)./2;
        % get bump full width at half maximum in pixel
        fwhm_pix=2.355*sigmaequiv;
        htot_pix=size(inputdata_temp,1); % screen height in pix
        htot_deg=89; % screen height in deg
        % get bump full width at half maximum in visual degrees
        fwhm_deg=htot_deg.*(fwhm_pix./htot_pix);
        % compute SF-equivalent
        SF_equiv=1/(2*fwhm_deg);
        % extract center
        X_center=fitParams_bestSTAlobe{target_area_idx}{neu_idx}(1);
        Y_center=fitParams_bestSTAlobe{target_area_idx}{neu_idx}(2);
        % store results
        SF_equiv_bestSTAlobe{target_area_idx}(neu_idx)=SF_equiv;
        X_center_bestSTAlobe{target_area_idx}(neu_idx)=X_center;
        Y_center_bestSTAlobe{target_area_idx}(neu_idx)=Y_center;
        % store classid
        classid_bestSTAlobe{target_area_idx}(neu_idx)=pop_classid{target_area_idx}(current_neuid);
        % visualize fit
        if boolplotfitdiagnostics
            % produce diagnostic figure
            fdiag=figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(2,1,1);
            imagesc(inputdata); colorbar; colormap(gray); caxis([-tempmaxval,tempmaxval]*0.9); axis equal
            hold on;
            numPoints=100;
            [ellipseX, ellipseY, ellipseX_center, ellipseY_center] = get_ellipse_from_2D_gaussian_fit(fitParams_bestSTAlobe{target_area_idx}{neu_idx},numPoints);
            plot(gca,ellipseX,ellipseY,'linewidth',2,'color',colortouse);
            scatter(gca,ellipseX_center,ellipseY_center,175,'*','Markerfacecolor',colortouse,'Markeredgecolor',colortouse)
            title('maxlobe STA frame')
            set(gca,'fontsize',12)
            subplot(2,1,2);
            imagesc(inputdata_temp); colorbar; colormap(gray); caxis([-tempmaxval,tempmaxval]*0.9); axis equal
            title('original STA frame')
            suptitle(['n good STA neu = ',num2str(neu_idx),' - ',areas{target_area_idx},' - SF equiv = ',num2str(round(SF_equiv,3))])
            set(gca,'fontsize',12)
            saveas(fdiag,[outfold,filesep,'RF_gaussianfit_diagnostics',filesep,areas{target_area_idx},'_RF_fit_diagnostics_n',num2str(neu_idx)],'jpg')
            print(fdiag,'-depsc','-painters',[outfold,filesep,'RF_gaussianfit_diagnostics',filesep,areas{target_area_idx},'_RF_fit_diagnostics_n',num2str(neu_idx),'.eps'])
            close all
        end
    end
    toc
end
% filter outliers
SF_equiv_bestSTAlobe_clean=SF_equiv_bestSTAlobe;
for iii=1:numel(SF_equiv_bestSTAlobe)
    SF_equiv_bestSTAlobe_clean{iii}=SF_equiv_bestSTAlobe_clean{iii}(SF_equiv_bestSTAlobe_clean{iii}<0.1);
end
% plot STA equiv. SF distributions  ------------------------------------
fighandplus=figure('units','normalized','outerposition',[0 0 1 1]);
distribtouse=SF_equiv_bestSTAlobe_clean;
inputpars.inputaxh=gca;
hold(inputpars.inputaxh,'on')
% set settings for violin distribution plotting
inputpars.boxplotwidth=0.4;%0.5;
inputpars.boxplotlinewidth=2;
inputpars.densityplotwidth=0.5;%0.5;
inputpars.yimtouse=[0,0.07];
% inputpars.yimtouse=[0,8];
inputpars.scatterjitter=inputpars.boxplotlinewidth*0.1;
inputpars.scatteralpha=0.15;
inputpars.scattersize=20;
inputpars.distralpha=0.5;
inputpars.xlabelstring=[];
inputpars.ylabelstring='STA equiv. SF';
inputpars.titlestring=[];
inputpars.boolscatteron=1;
inputpars.ks_bandwidth=0.0033;
inputpars.xlimtouse=[0,4]; %[-1,5];
% plot violins
inputadata.inputdistrs=distribtouse;
inputpars.n_distribs=numel(inputadata.inputdistrs);
inputpars.dirstrcenters=(1:inputpars.n_distribs);
inputpars.xtickslabelvector={'V1','LM','RL'};
inputpars.distrcolors{1}=([255,255,255]./255)*0.75;
inputpars.distrcolors{2}=([255,255,255]./255)*0.5;
inputpars.distrcolors{2}=([255,255,255]./255)*0.25;
inputaxh = plot_violinplot_PN_new(inputadata,inputpars);
hold on;
text(0.25,0.055,['avgV1 (#',num2str(numel(distribtouse{1})),') =',num2str(round(nanmean(distribtouse{1}),4))],'fontsize',12)
text(0.25,0.052,['avgLM (#',num2str(numel(distribtouse{2})),') =',num2str(round(nanmean(distribtouse{2}),4))],'fontsize',12)
text(0.25,0.049,['avgRL (#',num2str(numel(distribtouse{3})),') =',num2str(round(nanmean(distribtouse{3}),4))],'fontsize',12)
sgtitle(['automatic SF-equiv estimation (best frame) - all areas - all neurons CI>=',num2str(round(contrast_th_to_use,2))])
saveas(fighandplus,[outfold,filesep,'STA_best_frame_SF_eqiv_distr'],'jpg')
print(fighandplus,'-depsc','-painters',[outfold,filesep,'STA_best_frame_SF_eqiv_distr','.eps'])

%% analyze STA best frames RF centering fror pattern and component ----------------------

% set contrast threshold to use for fitting
contrast_th_to_use=7.5;
% set whether to plot fitting diagnostics
boolplotfitdiagnostics=0;
% set regularization to use for fitting
regpars=5*ones(1,5);
% loop over areas
for target_area_idx=1:numel(target_areas)
    tic
    % select good contrast STA neuron indeces
    goodcontrast_neuid{target_area_idx}=find(selectedframes_contrast{target_area_idx}>=contrast_th_to_use);
    % initialize current fit datastructures
    fitParams_bestSTAlobe{target_area_idx}=cell(1,numel(goodcontrast_neuid{target_area_idx}));
    fittedGaussian_bestSTAlobe{target_area_idx}=cell(1,numel(goodcontrast_neuid{target_area_idx}));
    SF_equiv_bestSTAlobe{target_area_idx}=NaN(1,numel(goodcontrast_neuid{target_area_idx}));
    X_center_bestSTAlobe{target_area_idx}=NaN(1,numel(goodcontrast_neuid{target_area_idx}));
    Y_center_bestSTAlobe{target_area_idx}=NaN(1,numel(goodcontrast_neuid{target_area_idx}));
    classid_bestSTAlobe{target_area_idx}=NaN(1,numel(goodcontrast_neuid{target_area_idx}));
    % loop over neurons with good STA in current area
    for neu_idx=1:numel(goodcontrast_neuid{target_area_idx})
        % get current neuron
        current_neuid=goodcontrast_neuid{target_area_idx}(neu_idx);
        % gest best STA frame for current neuron
        inputdata_temp=imresize(STA_bestframe{target_area_idx}(:,:,current_neuid),2);
        % take abs (keeping both excitatory and inhibitory lobe)
        inputdata_temp=abs(inputdata_temp);
        inputdata=inputdata_temp-nanmean(inputdata_temp);
        % light smoothing for improving fitting quality
        filter_size = [5 5];
        sigma = 1;
        filter = fspecial('gaussian', filter_size, sigma);
        inputdata = imfilter(inputdata, filter);
        % perform fitting
        [fitParams_bestSTAlobe{target_area_idx}{neu_idx},...
            fittedGaussian_bestSTAlobe{target_area_idx}{neu_idx}, ~]...
            = get_2D_gaussian_fit(inputdata,1.85*regpars);
        % extract center
        X_center=fitParams_bestSTAlobe{target_area_idx}{neu_idx}(1);
        Y_center=fitParams_bestSTAlobe{target_area_idx}{neu_idx}(2);
        % store results
        X_center_bestSTAlobe{target_area_idx}(neu_idx)=X_center;
        Y_center_bestSTAlobe{target_area_idx}(neu_idx)=Y_center;
        % store classid
        classid_bestSTAlobe{target_area_idx}(neu_idx)=pop_classid{target_area_idx}(current_neuid);
    end
    toc
end
% recover original neurn number
ns_per_class=cell(numel(target_areas),numel(cell_types_codes));
cell_types={'component','pattern','unclassified'};
for target_area_idx=1:length(target_areas)
    % set area
    are=target_areas{target_area_idx};
    % loop over original cell types
    for cell_type_idx=1:3
        cty=cell_types{cell_type_idx};
        rffit=load([outfold,filesep,are,'_',cty,'_prediction_explained_variance_addendum_2022.mat']);
        ns_per_class{target_area_idx,cell_type_idx}=rffit.ns_old;
    end
end
% set tolerance
toler=0;
% initialize out of screen boolean vector
bool_is_out_of_screen=cell(1,numel(target_areas));
% loop over areas
for target_area_idx=1:numel(target_areas)
    bool_is_out_of_screen{target_area_idx}=NaN(1,numel(X_center_bestSTAlobe{target_area_idx}));
    % loop over neurons
    for neu_idx=1:numel(X_center_bestSTAlobe{target_area_idx})
        % get FWHM ellipse
        parstouse=fitParams_bestSTAlobe{target_area_idx}{neu_idx}./2;
        numPoints=100;
        [ellipseX, ellipseY, ellipseX_center, ellipseY_center] = get_ellipse_from_2D_gaussian_fit(parstouse,numPoints);
        % determine if within or outside
        outX=logical(sum(not(and(ellipseX<size(pop_wSTA{target_area_idx},2)-0.5,ellipseX>0+0.5))));
        outY=logical(sum(not(and(ellipseY<size(pop_wSTA{target_area_idx},1)-0.5,ellipseY>0+0.5))));
        bool_is_out_of_screen{target_area_idx}(neu_idx)=or(outX,outY);
    end
end
% initialize out of screen fraction datastructures
frac_out_per_class=NaN(numel(target_areas),numel(cell_types_codes));
frac_out_per_class_num=NaN(numel(target_areas),numel(cell_types_codes));
frac_out_per_class_denom=NaN(numel(target_areas),numel(cell_types_codes));
STAs_per_class=cell(numel(target_areas),numel(cell_types_codes));
framecontrasts_per_class=cell(numel(target_areas),numel(cell_types_codes));
pop_rSTAs_per_class=cell(numel(target_areas),numel(cell_types_codes));
pop_wSTAs_per_class=cell(numel(target_areas),numel(cell_types_codes));
tot_frac_out_per_class=NaN(1,numel(cell_types_codes));
pop_TCs_per_class=cell(numel(target_areas),numel(cell_types_codes));
bool_is_out_of_screen_per_class=cell(numel(target_areas),numel(cell_types_codes));
% loop over areas
for target_area_idx=1:numel(target_areas)
    % loop over cell classes
    for cell_type_idx = 1:numel(cell_types_codes)
        % get currently selected neurons
        current_selected_id=find(classid_bestSTAlobe{target_area_idx}==cell_type_idx);
        % compute outside border fraction
        frac_out_per_class(target_area_idx,cell_type_idx)=sum(bool_is_out_of_screen{target_area_idx}(current_selected_id))./numel(bool_is_out_of_screen{target_area_idx}(current_selected_id));
        frac_out_per_class_num(target_area_idx,cell_type_idx)=sum(bool_is_out_of_screen{target_area_idx}(current_selected_id));
        frac_out_per_class_denom(target_area_idx,cell_type_idx)=numel(bool_is_out_of_screen{target_area_idx}(current_selected_id));
        STAs_per_class{target_area_idx,cell_type_idx}=pop_STAs{target_area_idx}(:,:,goodcontrast_neuid{target_area_idx}(current_selected_id));
        framecontrasts_per_class{target_area_idx,cell_type_idx}=selectedframes_contrast{target_area_idx}(goodcontrast_neuid{target_area_idx}(current_selected_id));
        pop_TCs_per_class{target_area_idx,cell_type_idx}=pop_TCs{target_area_idx}(:,goodcontrast_neuid{target_area_idx}(current_selected_id),:,:,:);
        % compute contrast weighted population STA per class
        weightingvid=permute(repmat(framecontrasts_per_class{target_area_idx,cell_type_idx},...
            [size(STAs_per_class{target_area_idx,cell_type_idx},1),1,size(STAs_per_class{target_area_idx,cell_type_idx},2)]),[1,3,2]);
        weightingvid=weightingvid./sum(framecontrasts_per_class{target_area_idx,cell_type_idx});
        curr_pop_wSTA=nansum(abs(STAs_per_class{target_area_idx,cell_type_idx}).*weightingvid,3);
        % compute unweighted population STA per class
        curr_pop_rSTA=nanmean(abs(STAs_per_class{target_area_idx,cell_type_idx}),3);
        % store population STA per class
        pop_rSTAs_per_class{target_area_idx,cell_type_idx}=curr_pop_rSTA;
        pop_wSTAs_per_class{target_area_idx,cell_type_idx}=curr_pop_wSTA;
        % index of neurons out of screen
        bool_is_out_of_screen_per_class{target_area_idx,cell_type_idx}=find(bool_is_out_of_screen{target_area_idx}(current_selected_id));
    end
end
% loop over cell classes
for cell_type_idx = 1:numel(cell_types_codes)
    % compute total fraction including all areas
    tot_frac_out_per_class(cell_type_idx)=sum((frac_out_per_class_num(:,cell_type_idx)))./sum((frac_out_per_class_denom(:,cell_type_idx)));
end

% % % % visualize fitteing procedure
% % % % set regularization to use for fitting
% % % regpars=5*ones(1,5);
% % % % loop over areas
% % % for target_area_idx=1:numel(target_areas)
% % %     % loop over cell classes
% % %     for cell_type_idx = 1:2%1:numel(cell_types_codes)
% % %         % loop over nezrons with good STA in current area
% % %         for neu_idx=1:size((STAs_per_class{target_area_idx,cell_type_idx}),3)
% % %             % gest best STA frame for current neuron
% % %             inputdata_orig=imresize(STAs_per_class{target_area_idx,cell_type_idx}(:,:,neu_idx),2);
% % %             % take abs
% % %             inputdata_temp=abs(inputdata_orig);
% % %             inputdata=inputdata_temp;
% % %             % light smoothing for improving fitting quality
% % %             filter_size = [15 15];
% % %             sigma = 1;
% % %             filter = fspecial('gaussian', filter_size, sigma);
% % %             inputdata = imfilter(inputdata, filter);
% % %             inputdata = inputdata./max(inputdata(:));
% % %             % perform fitting
% % %             [fitParams_per_class{target_area_idx}{neu_idx},...
% % %                 fittedGaussian_per_class{target_area_idx}{neu_idx}, ~] ...
% % %                 = get_2D_gaussian_fit(inputdata,1.85*regpars); %#ok<SAGROW>
% % %             % produce diagnostic figure
% % %             fdiag=figure('units','normalized','outerposition',[0 0 1 1]);
% % %             subplot(2,1,1);
% % %             maaax=quantile(inputdata(:),0.999);
% % %             miiin=quantile(inputdata(:),0.001);
% % %             [tempmaxval,tempid]=max([abs(miiin),abs(maaax)]);
% % %             imagesc(inputdata); colorbar; colormap(gray); caxis([0,tempmaxval]*1); axis equal
% % %             hold on;
% % %             numPoints=100;
% % %             [ellipseX, ellipseY, ellipseX_center, ellipseY_center] = get_ellipse_from_2D_gaussian_fit(fitParams_per_class{target_area_idx}{neu_idx},numPoints);
% % %             plot(gca,ellipseX,ellipseY,'linewidth',2,'color',colortouses{cell_type_idx});
% % %             scatter(gca,ellipseX_center,ellipseY_center,175,'*','Markerfacecolor',colortouses{cell_type_idx},'Markeredgecolor',colortouses{cell_type_idx})
% % %             title('maxlobe STA frame to be fitted')
% % %             set(gca,'fontsize',12)
% % %             subplot(2,1,2);
% % %             maaax=quantile(inputdata_orig(:),0.999);
% % %             miiin=quantile(inputdata_orig(:),0.001);
% % %
% % %             [tempmaxval,tempid]=max([abs(miiin),abs(maaax)]);
% % %             imagesc(inputdata_orig); hold on; colorbar; colormap(gray); caxis([-tempmaxval,tempmaxval]*1); axis equal
% % %             plot(gca,ellipseX,ellipseY,'linewidth',2,'color',colortouses{cell_type_idx});
% % %             scatter(gca,ellipseX_center,ellipseY_center,175,'*','Markerfacecolor',colortouses{cell_type_idx},'Markeredgecolor',colortouses{cell_type_idx})
% % %             title('original STA frame')
% % %             suptitle(['n good STA neu = ',num2str(neu_idx),' - ',areas{target_area_idx}])
% % %             pause(1)
% % %         end
% % %     end
% % % end

% plot average pop RF per class ------------------------------------
figuuuu=figure('units','normalized','outerposition',[0 0 1 1]);
labeltouse={'pattern','component'};
colortouses{1}=[50,200,0]./255;
colortouses{2}=[255,150,0]./255;
colortouses{3}=[0,0,0]./255;
% loop over areas
for target_area_idx=1:numel(target_areas)
    % loop over cell classes
    for cell_type_idx = 1:2
        curr_n_neu=size(STAs_per_class{target_area_idx,cell_type_idx},3);
        subplot(3,2,cell_type_idx+(target_area_idx-1)*2)
        imagesc(squeeze(pop_rSTAs_per_class{target_area_idx,cell_type_idx}));
        colorbar; colormap(gray);
        axis equal
        ylim([1-0.5,size(pop_rSTAs_per_class{target_area_idx,cell_type_idx},1)-0.5])
        xlim([1-0.5,size(pop_rSTAs_per_class{target_area_idx,cell_type_idx},2)-0.5])
        set(gca,'fontsize',12)
        title(['pop STA ',labeltouse{cell_type_idx},'( #',num2str(curr_n_neu),' )'],'color',colortouses{cell_type_idx})
    end
end
saveas(figuuuu,[outfold,filesep,'pop_STA_per_areas_and_class'],'jpg')
print(figuuuu,'-depsc','-painters',[outfold,filesep,'pop_STA_per_areas_and_class','.eps'])
% plot out of screen analysis results ------------------------------------
fighand00000000000=figure('units','normalized','outerposition',[0 0 1 1]);
colortouses{1}=[50,200,0]./255;
colortouses{2}=[255,150,0]./255;
colortouses{3}=[0,0,0]./255;
colortouse=([50,200,0]./255 + [255,150,0]./255)/2;
subplot(3,2,1)
imagesc(pop_wSTA{1}); colormap(gray); colorbar; axis square
hold on;
numPoints=100;
for neu_idx=1:numel(X_center_bestSTAlobe{1})
    if classid_bestSTAlobe{1}(neu_idx)~=3
        scatter(gca,X_center_bestSTAlobe{1}(neu_idx)./2,...
            Y_center_bestSTAlobe{1}(neu_idx)./2,200,'.',...
            'Markerfacecolor',colortouses{classid_bestSTAlobe{1}(neu_idx)},'Markeredgecolor',colortouses{classid_bestSTAlobe{1}(neu_idx)})
        parstouse=fitParams_bestSTAlobe{1}{neu_idx}./2;
        [ellipseX, ellipseY, ellipseX_center, ellipseY_center] = get_ellipse_from_2D_gaussian_fit(parstouse,numPoints);
        if bool_is_out_of_screen{1}(neu_idx)
            plot(gca,ellipseX,ellipseY,':','linewidth',3,'color',[colortouses{classid_bestSTAlobe{1}(neu_idx)},1]);
        else
            plot(gca,ellipseX,ellipseY,'linewidth',3,'color',[colortouses{classid_bestSTAlobe{1}(neu_idx)},1]);
        end
    end
end
title([target_areas{1},' grand avg (#comp=',...
    num2str(sum(classid_bestSTAlobe{1}==1)),...
    ' #patt=',num2str(sum(classid_bestSTAlobe{1}==2)),')']);
set(gca,'fontsize',12)
axis equal
xlim([1-0.5,size(pop_wSTA{1},2)+0.5])
ylim([1-0.5,size(pop_wSTA{1},1)+0.5])
subplot(3,2,3)
imagesc(pop_wSTA{2}); colormap(gray); colorbar; axis square
hold on;
numPoints=100;
for neu_idx=1:numel(X_center_bestSTAlobe{2})
    if classid_bestSTAlobe{2}(neu_idx)~=3
        scatter(gca,X_center_bestSTAlobe{2}(neu_idx)./2,...
            Y_center_bestSTAlobe{2}(neu_idx)./2,200,'.',...
            'Markerfacecolor',colortouses{classid_bestSTAlobe{2}(neu_idx)},'Markeredgecolor',colortouses{classid_bestSTAlobe{2}(neu_idx)})
        parstouse=fitParams_bestSTAlobe{2}{neu_idx}./2;
        [ellipseX, ellipseY, ellipseX_center, ellipseY_center] = get_ellipse_from_2D_gaussian_fit(parstouse,numPoints);
        if bool_is_out_of_screen{2}(neu_idx)
            plot(gca,ellipseX,ellipseY,':','linewidth',3,'color',[colortouses{classid_bestSTAlobe{2}(neu_idx)},1]);
        else
            plot(gca,ellipseX,ellipseY,'linewidth',3,'color',[colortouses{classid_bestSTAlobe{2}(neu_idx)},1]);
        end
    end
end
title([target_areas{2},' grand avg (#comp=',...
    num2str(sum(classid_bestSTAlobe{2}==1)),...
    ' #patt=',num2str(sum(classid_bestSTAlobe{2}==2)),')']);
set(gca,'fontsize',12)
axis equal
xlim([1-0.5,size(pop_wSTA{2},2)+0.5])
ylim([1-0.5,size(pop_wSTA{2},1)+0.5])
subplot(3,2,5)
imagesc(pop_wSTA{3}); colormap(gray); colorbar; axis square
hold on;
numPoints=100;
for neu_idx=1:numel(X_center_bestSTAlobe{3})
    if classid_bestSTAlobe{3}(neu_idx)~=3
        scatter(gca,X_center_bestSTAlobe{3}(neu_idx)./2,...
            Y_center_bestSTAlobe{3}(neu_idx)./2,200,'.',...
            'Markerfacecolor',colortouses{classid_bestSTAlobe{3}(neu_idx)},'Markeredgecolor',colortouses{classid_bestSTAlobe{3}(neu_idx)})
        parstouse=fitParams_bestSTAlobe{3}{neu_idx}./2;
        [ellipseX, ellipseY, ellipseX_center, ellipseY_center] = get_ellipse_from_2D_gaussian_fit(parstouse,numPoints);
        if bool_is_out_of_screen{3}(neu_idx)
            plot(gca,ellipseX,ellipseY,':','linewidth',3,'color',[colortouses{classid_bestSTAlobe{3}(neu_idx)},1]);
        else
            plot(gca,ellipseX,ellipseY,'linewidth',3,'color',[colortouses{classid_bestSTAlobe{3}(neu_idx)},1]);
        end
    end
end
title([target_areas{3},' grand avg (#comp=',...
    num2str(sum(classid_bestSTAlobe{3}==1)),...
    ' #patt=',num2str(sum(classid_bestSTAlobe{3}==2)),')']);
set(gca,'fontsize',12)
axis equal
xlim([1-0.5,size(pop_wSTA{3},2)+0.5])
ylim([1-0.5,size(pop_wSTA{3},1)+0.5])
subplot(3,2,[2,4,6])
barheightstouse=[frac_out_per_class(:,1),frac_out_per_class(:,2)];
colorstouse={colortouses{1},colortouses{1},colortouses{1},colortouses{2},colortouses{2},colortouses{2}};
barposition=[1,2,3,4+1,5+1,6+1];
hold on;
for jj=1:numel(barposition)
    bar(barposition(jj),1-barheightstouse(jj),'facecolor',colorstouse{jj},'edgecolor',colorstouse{jj})
end
ylabel('fraction')
xticks(barposition)
xticklabels({'comp. V1','comp. LM','comp. RL','patt. V1','patt. LM','patt. RL'})
title('fraction of cells within border')
text(0.25,1,['total faction comp. within =',num2str(round(1-tot_frac_out_per_class(1),2))],'fontsize',12)
text(0.25,0.97,['total faction patt. within =',num2str(round(1-tot_frac_out_per_class(2),2))],'fontsize',12)
ylim([0,1.1])
set(gca,'fontsize',12)
sgtitle(['population STA (abs) - all areas'])
saveas(fighand00000000000,[outfold,filesep,'pop_STA_per_areas_pc'],'jpg')
print(fighand00000000000,'-depsc','-painters',[outfold,filesep,'pop_STA_per_areas_pc','.eps'])