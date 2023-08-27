function [ F1z, F1F0 , plothandle, plotdata ] = get_modulation_index( psth, psth_duration, TF, boolplot )
% [ F1z, F1F0, plothandle, plotdata ] = get_modulation_index( psth, psth_duration, TF, boolplot )
% compute F1z and F1/F0 indeces given a psth

%% preprocess input ------------------------

% interpolate input psth
interpfactor=2;
n=length(psth);
v=psth;
x=0:psth_duration/n:(n-1)*psth_duration/n;
xq=0:psth_duration/(interpfactor*n):(n-1)*psth_duration/n;
psth=interp1(x,v,xq,'pchip');
N=length(psth);

%% compute F1z ------------------------

% demean signal
ds=psth;
baseline=mean(ds);
ds = ds-baseline;

% compute power spectrum
N_f=3*N;
T_s=psth_duration/N;
f_s=1/T_s;
maxlag=fix(N/3);
lagwindow='t';
[pow,f]=compute_pBT(ds, maxlag, N_f,f_s,lagwindow);

% find stimulus frequency vector index
fidx=find(f<=TF);
fidx=fidx(end);
spect=pow(fidx); sigspect=std(pow); meanspect=mean(pow);
F1z = (spect-meanspect)/sigspect;

%% compute F1/F0 ------------------------

% power based
ds=psth;
N_f=3*N;
T_s=psth_duration/N;
f_s=1/T_s;
maxlag=fix(N/3);
lagwindow='t';
[pow_bis,f]=compute_pBT(ds, maxlag, N_f,f_s,lagwindow);
fidx=find(f<=TF);
fidx=fidx(end);
F1F0 = pow_bis(fidx)/pow_bis(1);

%% store data for plotting -------------

plotdata.N_f=N_f;
plotdata.f=f;
plotdata.fidx=fidx;
plotdata.pow=pow;
plotdata.meanspect=meanspect;
plotdata.sigspect=sigspect;
plotdata.TF=TF;
plotdata.F1F0=F1F0;
plotdata.F1z=F1z;

%% plot results ------------------------

if boolplot
    
    plothandle=figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1);
    psth_x=0:T_s:(length(psth)-1)*T_s;
    plot(psth_x,psth,'k','LineWidth',3)
    hlbx=get(gca,'Xlabel');
    set(hlbx,'String','time [s]','FontWeight','bold','FontSize',10,'color','k')
    hlby=get(gca,'Ylabel');
    set(hlby,'String','unit activity','FontWeight','bold','FontSize',10,'color','k')
    title(['Unit response over time']);
    axis square
    set(gca,'fontsize',12);
    subplot(1,2,2);
    plot(f(1:round(N_f/2)+1),pow(1:round(N_f/2)+1),'k','LineWidth',3);
    hold on
    plot(f(fidx),pow(fidx),'ob','LineWidth',6);
    plot(f(1:round(N_f/2)+1),(meanspect+sigspect)*ones(size(pow(1:round(N_f/2)+1))),'--b','LineWidth',1);
    plot(f(1:round(N_f/2)+1),(meanspect-sigspect)*ones(size(pow(1:round(N_f/2)+1))),'--b','LineWidth',1);
    plot(f(1:round(N_f/2)+1),(meanspect)*ones(size(pow(1:round(N_f/2)+1))),'.-b','LineWidth',0.5);
    ylimit=get(gca,'ylim');
    xlimit=get(gca,'xlim');
    ttt=text(0.6*xlimit(2),0.85*ylimit(2),['target TF = ',num2str(TF),' Hz'],'FontSize',14); %#ok<NASGU>
    set(gca,'FontSize',10);
    hlabelx=get(gca,'Xlabel');
    set(hlabelx,'String','f [Hz]','FontWeight','bold','FontSize',10,'color','k')
    hlabely=get(gca,'Ylabel');
    set(hlabely,'String','PSD','FontWeight','bold','FontSize',10,'color','k')
    title(['Power spectrum (F1z = ',num2str(F1z),', F1F0 = ',num2str(F1F0),')']);  
    hold off
    axis square
    set(gca,'fontsize',12);
else
    plothandle=[];
end

end

