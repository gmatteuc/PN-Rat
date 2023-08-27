function pred_err = prediction_error(input_PSTHs,input_S,input_dirs,input_pdir,...
    alpha,beta,input_F,sr,countingwindowlim)

% get prediction
[obs_fr,pred_fr,obs_count,pred_count,~] =...
    predict_neuron_responses(input_PSTHs,input_S,input_dirs,...
    alpha,beta,input_F,sr,countingwindowlim);

% normalize predicted and observed data to 1
obs_fr=obs_fr./max(obs_fr,[],1);
pred_fr=pred_fr./max(pred_fr,[],1);
obs_count=obs_count./max(obs_count);
pred_count=pred_count./max(pred_count);

% get goodness of fit
[explvar_tuning,~,~,~] = get_neuron_responses_prediction_gof(obs_fr,obs_count,pred_fr,pred_count,input_dirs==input_pdir);
% get overall prediction error
unexplvar_tuning=1-explvar_tuning;
pred_err=unexplvar_tuning;%(rmse_tuning+rmse_dynamics)/2;

end

