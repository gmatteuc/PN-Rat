function [resp_dy,resp] = blobby_filtering_PN(filter,stims)
resp_dy=NaN(size(stims{1},3),numel(stims));
resp=NaN(1,numel(stims));
for i=1:numel(stims)
    for j=1:size(stims{i},3)
        curr_frame=stims{i}(:,:,j);
        curr_product=filter(:)'*curr_frame(:);
        resp_dy(j,i)=curr_product;
    end
    resp_dy(:,i)=max(resp_dy(:,i),zeros(size(resp_dy(:,i))));
    resp(i)=max(resp_dy(:,i));
end
end

