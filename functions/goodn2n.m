function [ n ] = goodn2n( area,goodn )

% load tuning analysis results and indexing
load(['D:\Backups\Personal_bk\PN_acute_analysis\processed_data',filesep,'Indexing.mat'])
% load neuron selection results
load(['D:\Backups\Personal_bk\PN_acute_analysis\processed_data',filesep,area,'_selected_neurons',filesep,'selected_population_',area,'.mat'])
n=selectedsi{goodn}(1);

end

