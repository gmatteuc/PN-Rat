function [ PI, PIgen, Zp, Zc, Rp, Rc, Zpgen, Zcgen, Rpgen, Rcgen ] = get_plaid_analysis_PN( nn, sf, tf, tuning_curve)

%[ PI, PIgen, Zp, Zc, Rp, Rc, Zpgen, Zcgen, Rpgen, Rcgen ] = get_plaid_analysis_PN( nn, sf, tf, tuning_curve)
%
%
%Perform plaid analysis

%--------------------------------------------------------------------------

% set some pars
pars = set_pars_PN();
SF=pars.stimPars.SF;
TF=pars.stimPars.TF;

% slice for selected bitcodes - gratings
tuning_curve_gr=tuning_curve(:,nn,SF==sf,TF==tf,1);

% slice for selected bitcodes - plaids
tuning_curve_pl=tuning_curve(:,nn,SF==sf,TF==tf,2);

% compute patterna and component predictions
[ PI, PIgen, Zp, Zc, Rp, Rc, Zpgen, Zcgen, Rpgen, Rcgen ]  = get_pattern_index( tuning_curve_gr,tuning_curve_pl );

end