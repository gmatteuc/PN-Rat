function [obs_fr,pred_fr,obs_count,pred_count,fr_time] = predict_neuron_responses(input_psths,input_stimuli,input_dir,alpha,beta,gamma,sr,counting_window_lim)

counting_window_span=counting_window_lim(2)-counting_window_lim(1);
% initialize storage variables
temp_obs_COUNT=NaN(1,length(input_dir));
temp_pred_COUNT=NaN(1,length(input_dir));
temp_obs_FR=NaN(sr.*(counting_window_span),length(input_dir));
temp_pred_FR=NaN(sr.*(counting_window_span),length(input_dir));
temp_pred_FR_time=cell(1,length(input_dir));

for dir_idx=1:length(input_dir)
    
    % get observed psth
    current_psth = input_psths{dir_idx};
    
    % smooth observed psth to match timescale of prediction
    gaussFilter = gausswin(10,1/0.5);
    gaussFilter = gaussFilter / sum(gaussFilter);
    current_psth=conv(current_psth, gaussFilter,'same');
    % resample observed psth
    current_psth=squeeze(current_psth);
    new_samplenum=sr.*(counting_window_span);
    old_samplenum=length(current_psth);
    current_psth_resampled = resample(current_psth,new_samplenum,old_samplenum);
    temp_obs_FR(:,dir_idx)=current_psth_resampled;
    
    % get observed count
    temp_obs_COUNT(dir_idx)=nansum(current_psth_resampled);
    
    % get current direction stimulus
    current_stimulus = input_stimuli{dir_idx};
    
    % get predicted psth
    input_filter=gamma;
    [current_pred_fr,temp_pred_COUNT(dir_idx),temp_pred_FR_time{dir_idx}]=...
        neural_filtering(input_filter,current_stimulus,alpha,beta,sr,counting_window_lim);
    % cut predicted psth
    tinf=find(temp_pred_FR_time{dir_idx}<1);
    tsup=find(temp_pred_FR_time{dir_idx}>0);
    tind=intersect(tinf,tsup);
    temp_pred_FR(:,dir_idx)=current_pred_fr(tind);
    
end

% get back prediction time
FR_time=temp_pred_FR_time{1};
tinf=find(FR_time<1);
tsup=find(FR_time>0);
tind=intersect(tinf,tsup);
fr_time=FR_time(tind);

% assign output variables
obs_fr=temp_obs_FR;
pred_fr=temp_pred_FR;
obs_count=temp_obs_COUNT;
pred_count=temp_pred_COUNT;

end

