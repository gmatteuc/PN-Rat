close all
clear all
clc

% README: this test function creates a dummy psth with a periodic 
% (i.e. simple-cell-like) and a constant (i.e. complex-cell-like) 
% component and calls the function "get_modulation_index" on it 
% to compute F1z and F1F0 with Blackman-Tukey spectral estimation

% generate test psth
psth=zeros(1,46);
psth(10:end)=3; % constant component 
psth(10:5:45)=15; % periodic component 
psth=smooth(psth,'sgolay');
psth(psth<0)=0;
psth_dur=1.5; % s
TF_target=6; % Hz

% get modulation indeces
[ F1z, F1F0 , plothandle ] = get_modulation_index( psth, psth_dur, TF_target, 1 );