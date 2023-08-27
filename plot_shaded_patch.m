function plot_shaded_patch(ax,bool_input,patch_limits,colortouse,alphatouse)

% get simply connected regions from boolean
region_starts=find(diff([0,bool_input,0])==1);
region_stops=find(diff([0,bool_input,0])==-1);
region_starts=region_starts-1;
region_starts(region_starts==0)=1;
region_stops=region_stops-1;

% loop over regions
for i=1:numel(region_starts)
    
    hold(ax,'on');
    
    if not(isempty(region_starts(i):region_stops(i)))
        
        % get current region boolean input
        current_bool_input=zeros(size(bool_input));
        current_bool_input(region_starts(i):region_stops(i))=1;
        
        % draw patch
        x = find(current_bool_input);
        uE = patch_limits(2).*ones(size(x));
        lE = patch_limits(1).*ones(size(x));
        yP = [lE,fliplr(uE)];
        xP = [x,fliplr(x)];
        idx2rm=isnan(yP);
        xP(idx2rm)=[];
        yP(idx2rm)=[];
        patch(ax,xP,yP,1,'facecolor',colortouse,...
            'edgecolor','none',...
            'facealpha',alphatouse);
        
    end
    
end

end

