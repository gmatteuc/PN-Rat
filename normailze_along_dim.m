function output_data = normailze_along_dim(input_data,input_d)

% get the number of elements along target dimension
num_target_elements = size(input_data, input_d);
% initialize an array to store the normalized input_data
output_data = zeros(size(input_data));
% get all indexes for normalization
input_data_all_idx=cell(1,numel(size(input_data)));
for input_data_dim_idx=1:numel(size(input_data))
    input_data_all_idx{input_data_dim_idx}=':';
end
% loop through elements along target dimension and normalize
for target_element_idx = 1:num_target_elements
    % get current indexes for normalization
    input_data_current_idx=input_data_all_idx;
    input_data_current_idx{input_d}=target_element_idx;
    % get current values
    current_vals=input_data(input_data_current_idx{:});
    % compute current max
    current_max=max(current_vals(:));
    % assigne normalized values to output
    output_data(input_data_current_idx{:})=current_vals./current_max;
end

end