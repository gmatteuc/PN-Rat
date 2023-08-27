function [] = plot_rasters_good( selectedsi, gooddir, gooddir_p, goomlabel, goodprds)

% plot_rasters( selectedsi, gooddir )
%
% plot grating and plaid rasters
%-----------------------------------------------------------------------

% set pars
pars=set_pars_PN;
data_folder=pars.processed_data_folder;
load(fullfile(data_folder,'Tuning.mat'));
load(fullfile(data_folder,'Indexing.mat'));
SF=pars.stimPars.SF;
TF=pars.stimPars.TF;
DIR=pars.stimPars.DIR;
stimulustypes={'grating','plaid'};
pre_time=pars.stimPars.pre_time;
post_time=pars.stimPars.post_time;
stim_time=pars.stimPars.stim_time;
post_delay=pars.stimPars.post_delay;
ntrials=pars.stimPars.ntrials;

cell_type='';
istoplot=0;

for goondn=1:size(selectedsi,1)
    
    if (goomlabel(goondn)*goodprds(goondn))==1
        istoplot=1;
        cell_type='component';
    elseif (goomlabel(goondn)*goodprds(goondn))==2
        istoplot=1;
        cell_type='pattern';
    else
        istoplot=0;
        cell_type='unclassified';
    end
    
    
    if istoplot
    
    for stimidx=1:2
        
        % switch colors
        if stimidx==1
            colo=[0,0,0.9]; % grating
        else
            colo=[0,0.7,0]; % plaids
        end
        
        % get current neuron parameters
        stimulustype=stimulustypes{stimidx};
        nn=selectedsi{goondn}(1);
        SFind=selectedsi{goondn}(2);
        sf=SF(SFind);
        TFind=selectedsi{goondn}(3);
        tf=TF(TFind);
        % switch colors
        if stimidx==1
            dir=DIR(gooddir(goondn)); % grating
        else
            dir=DIR(gooddir_p(goondn)); % plaids
        end
        
        % retrive within session index
        [intn,sessionname]=n2intn(nn);
        % get raster data
        finame=['SPIKEMAT_',sessionname];
        S=load(finame);
        
            % set plotname
    fname=['raster mosaic ',cell_type,' ',stimulustype,' SF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),' n',num2str(nn),''];
    fname=strrep(fname,'.','');
    fname=strrep(fname,' ','_');
    
    % initialize figure
    f1 = figure;
    set(f1,'Position',[10,10,1500,1000]);
    poscount1=0;
    poscount2=0;
    
    for dind=1:length(DIR)
        
        % select direction
        dir=DIR(dind);
        % select bitcode
        selected_idx = get_indexes_PN( sf, tf, dir, stimulustype );
        
        % initialize subplot
        sb1=subplot(666,666,666);
        hold on
        
        spp=cell(1,ntrials);
        for trial=1:ntrials
            if trial<=size(S.SPIKEtimes_bis( selected_idx, intn, :),3)
                sp=S.SPIKEtimes_bis{ selected_idx, intn, trial}-pre_time;
            else
                sp=[];
            end
            if not(isempty(sp))
                plot(sp,trial,'.k', 'MarkerSize',20)
            else
            end
            box on
            spp{trial}=sp;
        end
        xlim([-pre_time,post_time+stim_time]);
        ylim([-5,ntrials+5]);
        plot([0,0],[-5,ntrials+5],'--', 'LineWidth',3,'Color',colo)
        plot([0.9,0.9],[-5,25],'--', 'LineWidth',3,'Color',colo)
        text(-0.5,22,['DIR = ',num2str(dir),' ',' number of spikes = ',num2str(round(ntrials*response_fr(dir==DIR,nn,sf==SF,tf==TF,stimidx))),...
            ' (',num2str(round(ntrials*spontaneous_fr(dir==DIR,nn,sf==SF,tf==TF,stimidx))),')'],'FontSize',12);
        set(gca,'FontSize',10);
        hlabelx=get(gca,'Xlabel');
        set(hlabelx,'String','time from onset [s]','FontWeight','bold','FontSize',12,'color','k')
        hlabely=get(gca,'Ylabel');
        set(hlabely,'String','trial number','FontWeight','bold','FontSize',12,'color','k')
            
        % displace subplot to correct position
        hold off
        set(sb1,'Position',[.00+0.23*(poscount1)+0.07,.01+0.32*(poscount2)+0.1,0.20,0.20]);
        poscount1=poscount1+1;
        if poscount1==4 % change row after 4 subplots
            poscount2=poscount2+1;
            poscount1=0;
        end
        
    end
    
     % ---------------------------------------------------------------------
     
     % save plot
     saveas(f1,fname,'jpg')
     close all
    end
    
    message=['\nNeuron ',num2str(intn),' ',cell_type,' raster mosaic produced\n'];
    fprintf(message)
    clear S
    
    else
        
    message=['\n... plot skipped ...\n'];
    fprintf(message)
    
    end
    
end
end

