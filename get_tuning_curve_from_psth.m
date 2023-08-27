function [obs_fr,obs_count] = get_tuning_curve_from_psth(input_psths,sr,counting_window_lim)

counting_window_span=counting_window_lim(2)-counting_window_lim(1);
% initialize storage variables
temp_obs_COUNT=NaN(1,length(input_dir));
temp_obs_FR=NaN(sr.*(counting_window_span),length(input_dir));

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
    
end

% assign output variables
obs_fr=temp_obs_FR;
obs_count=temp_obs_COUNT;

end

