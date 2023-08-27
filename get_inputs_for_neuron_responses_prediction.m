function [input_PSTHs, input_RASTERs, input_STIMs, input_F, input_S] = get_inputs_for_neuron_responses_prediction(...
    input_interp_factor,input_neuron_num,input_sf,input_tf,input_dirs,input_stimtype,input_r2_per_frame,input_STA,r2_th)

% initialize current filter
input_filter = zeros(...
    size(input_STA,1)*input_interp_factor,...
    size(input_STA,2)*input_interp_factor,...
    size(input_STA,3));
% loop over frames of current filter
for ff=1:size(input_filter,3)
    % get current STA frame
    input_fr = input_STA(:,:,ff);
    input_fr_interp = interpolate_RF_frame( input_fr, input_interp_factor );
    % check if current frame is good
    is_frame_good_bool=input_r2_per_frame(ff)<=r2_th;
    if is_frame_good_bool
        % stack into filter matrix
        input_filter(:,:,ff) = zeros(size(input_filter,1),size(input_filter,2));
    else
        % stack into filter matrix
        input_filter(:,:,ff) = input_fr_interp;
    end
end

% get observed psths and visual stimuli
input_PSTHs=cell(1,length(input_dirs));
input_RASTERs=cell(1,length(input_dirs));
input_STIMs=cell(1,length(input_dirs));
for dir_idx=1:length(input_dirs)
    % get current direction
    current_direction=input_dirs(dir_idx);
    % get observed psth
    [input_PSTHs{dir_idx},input_RASTERs{dir_idx},~] = get_psth_PN(input_neuron_num,input_sf,input_tf,current_direction,input_stimtype);
    % get current direction stimulus
    input_STIMs{dir_idx} = get_stimulus_PN(input_sf,input_tf,current_direction,input_stimtype);
end

% reshaping filter for fast convolution
STAf=input_filter-mean(mean(mean(input_filter)));
STAf=STAf./(max(STAf(:))-min(STAf(:)));
flipped_STAf=flip(STAf,3);
% unroll filter
input_F=flipped_STAf(:);

% reshaping stimulus for fast convolution
num_conv_steps=size(input_STIMs{1},3)-size(input_filter,3);
num_stimblock_elements=numel(input_F);
input_S=cell(1,length(input_dirs));
for dir_idx=1:length(input_dirs)
    input_S{dir_idx}=NaN(num_stimblock_elements,num_conv_steps);
    for kk=1:size(input_STIMs{dir_idx},3)-size(STAf,3)
        % unroll current stimblock
        input_S{dir_idx}(:,kk)=single(reshape(input_STIMs{dir_idx}(:,:,(1:size(STAf,3))+(kk-1)),[1,num_stimblock_elements]));
    end
end

end

