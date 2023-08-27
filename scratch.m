

        % TODO: improve prediction  working on neural filtering function,
        % why peks offset? oprimize nonlinearity for each neuron? ensure
        % orientation continuity og gabor? ... neural filtering function
        % wrongly reverses the time? why desite sharpening zc didn-t
        % improve much?
        
% add info on gabor goodness orientation continuity og gabor?
        % already think to producing the data for the decoding: you will
        % need psth for all direction (pref stim) and original one to
        % estimate baseline noise level by having a costant background rate
        % obtained by averaging the residual psth (average across dirs and
        % rectify)

                
toc
                
                pred_FR_time=cell(1,length(DIR));
                parfor dir_idx=1:length(DIR)
                    
                    % get current direction
                    current_direction=DIR(dir_idx);
                    
                    % get observed psth
                    [current_psth,~,~] = get_psth_PN(n,pSF,pTF,current_direction,stimulustype);
                    % smooth observed psth to match timescale of prediction
                    gaussFilter = gausswin(10,1/0.5);
                    gaussFilter = gaussFilter / sum(gaussFilter);
                    current_psth=conv(current_psth, gaussFilter,'same');
                    % resample observed psth
                    current_psth=squeeze(current_psth);
                    new_samplenum=length(tind);
                    old_samplenum=length(current_psth);
                    current_psth_resampled = resample(current_psth,new_samplenum,old_samplenum);
                    obs_FR(:,dir_idx,i,k)=current_psth_resampled;
                    
                    % get current direction stimulus
                    current_stimulus = get_stimulus_PN(pSF,pTF,current_direction,stimulustype);
                    
                    % get predicted psth
                    [current_pred_fr,pred_COUNT(dir_idx,i,k),pred_FR_time{dir_idx}]=...
                        neural_filtering(current_filter,current_stimulus,alpha,beta,sr,countingwindowlim);
                    % cut predicted psth
                    tinf=find(pred_FR_time{dir_idx}<1);
                    tsup=find(pred_FR_time{dir_idx}>0);
                    tind=intersect(tinf,tsup);
                    pred_FR(:,dir_idx,i,k)=current_pred_fr(tind);
                    
                end

                % get back prediction time
                FR_time=pred_FR_time{1};
                tinf=find(FR_time<1);
                tsup=find(FR_time>0);
                tind=intersect(tinf,tsup);
                FR_time=FR_time(tind);
                % get observed tuning curve
                obs_COUNT(:,i,k)=tuning_curve(:,n,pSF==SF,pTF==TF,k);
                
                % output message
                fprintf(['filter response at ',stimulustype,' ( all DIRs',' SF=',num2str(pSF),' TF=',num2str(pTF),' ) neuron ',num2str(n),' computed ...\n']);
                toc