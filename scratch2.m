

% % needed to recompute pattern index
% input_PSTHs= get_psth_PN(input_neuron_num,input_sf,input_tf,current_direction,input_stimtype);
% sr=videoframerate;
% countingwindowlim=[0,1];
% [obs_fr,obs_count] = get_tuning_curve_from_psth(input_psths,sr,counting_window_lim)
% % get tuning curves
% tuning_curve_grating_P=pred_COUNT(:,i,1);
% tuning_curve_plaid_P=pred_COUNT(:,i,2);
% % perform partial correlation analysis to get predicted index values
% [ PI_P(i), ~, Zp_P(i), Zc_P(i), Rp_P(i), Rc_P(i), ~, ~, ~, ~ ] =...
%     get_pattern_index( tuning_curve_grating_P,tuning_curve_plaid_P );