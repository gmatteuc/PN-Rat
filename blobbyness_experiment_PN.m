close all
clear all
clc

%% get Grating and Plaid stimuli to use for blobbyness experiment

% set number of points
npts=500;
% set grating path
stims_path='D:\Backups\Personal_bk\PN_acute_analysis\stims';
% sf to use
sftouse='0.02';
% get names of folders containing stimuli
switch sftouse
    case '0.04'
        % set grating names
        grating_foldnames{1}='0136_grating_0_0.04_6';
        grating_foldnames{2}='0137_grating_30_0.04_6';
        grating_foldnames{3}='0138_grating_60_0.04_6';
        grating_foldnames{4}='0139_grating_90_0.04_6';
        grating_foldnames{5}='0140_grating_120_0.04_6';
        grating_foldnames{6}='0141_grating_150_0.04_6';
        grating_foldnames{7}='0142_grating_180_0.04_6';
        grating_foldnames{8}='0143_grating_210_0.04_6';
        grating_foldnames{9}='0144_grating_240_0.04_6';
        grating_foldnames{10}='0145_grating_270_0.04_6';
        grating_foldnames{11}='0146_grating_300_0.04_6';
        grating_foldnames{12}='0147_grating_330_0.04_6';
        % set plaid names
        plaid_foldnames{1}='0386_plaid_0_0.04_6';
        plaid_foldnames{2}='0387_plaid_30_0.04_6';
        plaid_foldnames{3}='0388_plaid_60_0.04_6';
        plaid_foldnames{4}='0389_plaid_90_0.04_6';
        plaid_foldnames{5}='0390_plaid_120_0.04_6';
        plaid_foldnames{6}='0391_plaid_150_0.04_6';
        plaid_foldnames{7}='0392_plaid_180_0.04_6';
        plaid_foldnames{8}='0393_plaid_210_0.04_6';
        plaid_foldnames{9}='0394_plaid_240_0.04_6';
        plaid_foldnames{10}='0395_plaid_270_0.04_6';
        plaid_foldnames{11}='0396_plaid_300_0.04_6';
        plaid_foldnames{12}='0397_plaid_330_0.04_6';
    case '0.02'
        % set grating names
        grating_foldnames{1}='0124_grating_0_0.02_6';
        grating_foldnames{2}='0125_grating_30_0.02_6';
        grating_foldnames{3}='0126_grating_60_0.02_6';
        grating_foldnames{4}='0127_grating_90_0.02_6';
        grating_foldnames{5}='0128_grating_120_0.02_6';
        grating_foldnames{6}='0129_grating_150_0.02_6';
        grating_foldnames{7}='0130_grating_180_0.02_6';
        grating_foldnames{8}='0131_grating_210_0.02_6';
        grating_foldnames{9}='0132_grating_240_0.02_6';
        grating_foldnames{10}='0133_grating_270_0.02_6';
        grating_foldnames{11}='0134_grating_300_0.02_6';
        grating_foldnames{12}='0135_grating_330_0.02_6';
        % set plaid names
        plaid_foldnames{1}='0374_plaid_0_0.02_6';
        plaid_foldnames{2}='0375_plaid_30_0.02_6';
        plaid_foldnames{3}='0376_plaid_60_0.02_6';
        plaid_foldnames{4}='0377_plaid_90_0.02_6';
        plaid_foldnames{5}='0378_plaid_120_0.02_6';
        plaid_foldnames{6}='0379_plaid_150_0.02_6';
        plaid_foldnames{7}='0380_plaid_180_0.02_6';
        plaid_foldnames{8}='0381_plaid_210_0.02_6';
        plaid_foldnames{9}='0382_plaid_240_0.02_6';
        plaid_foldnames{10}='0383_plaid_270_0.02_6';
        plaid_foldnames{11}='0384_plaid_300_0.02_6';
        plaid_foldnames{12}='0385_plaid_330_0.02_6';
end

% load gratings
gratings=cell(1,numel(grating_foldnames));
for i=1:numel(grating_foldnames)
    tic
    current_path=[stims_path,filesep,grating_foldnames{i}];
    current_dir=dir([current_path,filesep,'*.bmp*']);
    gratings{i}=NaN(npts,npts,numel(current_dir));
    for j=1:numel(current_dir)
        curr_img_path=[current_dir(j).folder,filesep,current_dir(j).name];
        curr_img=imresize(2*((double(imread(curr_img_path))./255)-0.5),[npts,npts]);
        gratings{i}(:,:,j)=curr_img;
    end
    toc
end
% load plaids
plaids=cell(1,numel(plaid_foldnames));
for i=1:numel(plaid_foldnames)
    tic
    current_path=[stims_path,filesep,plaid_foldnames{i}];
    current_dir=dir([current_path,filesep,'*.bmp*']);
    plaids{i}=NaN(npts,npts,numel(current_dir));
    for j=1:numel(current_dir)
        curr_img_path=[current_dir(j).folder,filesep,current_dir(j).name];
        curr_img=imresize(2*((double(imread(curr_img_path))./255)-0.5),[npts,npts]);
        plaids{i}(:,:,j)=curr_img;
    end
    toc
end

%% define filters to use for blobbyness experiment

% define fixed parameters for the Gabor filter
lambda = 50;
theta = deg2rad(90);
sigma = 0.2*lambda;
delta = deg2rad(0);

% create  meshgrid for the Gabor filters
[X, Y] = meshgrid(linspace(-50,50,npts), linspace(-50,50,npts));
% create derived meshgrid rotated of desired angle
Xprime=X.*cos(delta)+Y.*sin(delta);
Yprime=-X.*sin(delta)+Y.*cos(delta);

% evaluate Gabor function on meshgrid to produece filters with different aspect ratio
gamma = 0.4;
aspect1 = gamma;
Gabor1 = exp(-((Xprime.^2 + gamma^2 * Yprime.^2)/(2 * sigma^2))) .* cos(2*pi*Xprime/lambda + theta);
gamma = 0.60;
aspect2 = gamma;
Gabor2 = exp(-((Xprime.^2 + gamma^2 * Yprime.^2)/(2 * sigma^2))) .* cos(2*pi*Xprime/lambda + theta);
gamma = 0.7;
aspect3 = gamma;
Gabor3 = exp(-((Xprime.^2 + gamma^2 * Yprime.^2)/(2 * sigma^2))) .* cos(2*pi*Xprime/lambda + theta);

%% run blobbyness experiment (compute filters responses)

% set stim directions
stims_dir=0:30:330;

% colors
col1=[50,200,0]./255;
col2=[255,150,0]./255;

% get filter 1 gratings responses
stims=gratings;
filter=Gabor1-nanmean(Gabor1(:));
[gratings_resp_dy_f1,gratings_resp_f1] = blobby_filtering_PN(filter,stims);
% get filter 1 plaids responses
stims=plaids;
filter=Gabor1-nanmean(Gabor1(:));
[plaids_resp_dy_f1,plaids_resp_f1] = blobby_filtering_PN(filter,stims);

% get filter 2 gratings responses
stims=gratings;
filter=Gabor2-nanmean(Gabor2(:));
[gratings_resp_dy_f2,gratings_resp_f2] = blobby_filtering_PN(filter,stims);
% get filter 2 plaids responses
stims=plaids;
filter=Gabor2-nanmean(Gabor2(:));
[plaids_resp_dy_f2,plaids_resp_f2] = blobby_filtering_PN(filter,stims);

% get filter 3 gratings responses
stims=gratings;
filter=Gabor3-nanmean(Gabor3(:));
[gratings_resp_dy_f3,gratings_resp_f3] = blobby_filtering_PN(filter,stims);
% get filter 3 plaids responses
stims=plaids;
filter=Gabor3-nanmean(Gabor3(:));
[plaids_resp_dy_f3,plaids_resp_f3] = blobby_filtering_PN(filter,stims);

%% visualize blobbyness experiment results (visualize filters responses)

% plot blobbyness effect visualization
f0=figure('units','normalized','outerposition',[0 0 1 1]);
interpfact=50;
% aspect ratio 1 ----------------------------------------------
subplot(3,4,1)
imagesc(Gabor1); colormap('gray'); colorbar;
caxis([-max(abs(Gabor1(:))),+max(abs(Gabor1(:)))]);
axis square;
xlabel('X');
ylabel('Y');
title(['Gabor filter ( aspect ratio = ',num2str(aspect1),' )']);
set(gca,'fontsize',12)
subplot(3,4,2)
imagesc(gratings_resp_dy_f1); colormap(gray); colorbar;
ylabel('timebin #')
xlabel('direction (?)')
title('Grating')
xticks(1:numel(stims_dir)); xtickangle(45);
xticklabels(split(cellstr(num2str(stims_dir))))
axis square
set(gca,'fontsize',12)
subplot(3,4,3)
imagesc(plaids_resp_dy_f1); colormap(gray); colorbar;
ylabel('timebin #')
xlabel('direction (?)')
title('Plaid')
xticks(1:numel(stims_dir)); xtickangle(45);
xticklabels(split(cellstr(num2str(stims_dir))))
axis square
set(gca,'fontsize',12)
subplot(3,4,4)
hold on;
int_dirs=linspace(min(stims_dir),max(stims_dir),interpfact);
int_grating=interp1(stims_dir,gratings_resp_f1,int_dirs,'spline');
int_plaid=interp1(stims_dir,plaids_resp_f1,int_dirs,'spline');
int_grating=int_grating./max(int_grating(:));
int_plaid=int_plaid./max(int_plaid(:));
plot(int_dirs,int_grating,'color',col1,'linewidth',3)
plot(int_dirs,int_plaid,'color',col1/2,'linewidth',3)
plot([180,180],get(gca,'ylim'),':','color',[0.5,0.5,0.5],'linewidth',2)
xlim([180-90,180+90])
title('Grating and Plaid tuning curves')
ylim([0,1.1]);
ylabel('max response (normalized)');
xlabel('direction (?)')
set(gca,'fontsize',12)
% aspect ratio 2 ----------------------------------------------
subplot(3,4,5)
imagesc(Gabor2); colormap('gray'); colorbar;
caxis([-max(abs(Gabor2(:))),+max(abs(Gabor2(:)))]);
axis square;
xlabel('X');
ylabel('Y');
title(['Gabor filter ( aspect ratio = ',num2str(aspect2),' )']);
set(gca,'fontsize',12)
subplot(3,4,6)
imagesc(gratings_resp_dy_f2); colormap(gray); colorbar;
ylabel('timebin #')
xlabel('direction (?)')
title('Grating')
xticks(1:numel(stims_dir)); xtickangle(45);
xticklabels(split(cellstr(num2str(stims_dir))))
axis square
set(gca,'fontsize',12)
subplot(3,4,7)
imagesc(plaids_resp_dy_f2); colormap(gray); colorbar;
ylabel('timebin #')
xlabel('direction (?)')
title('Plaid')
xticks(1:numel(stims_dir)); xtickangle(45);
xticklabels(split(cellstr(num2str(stims_dir))))
axis square
set(gca,'fontsize',12)
subplot(3,4,8)
hold on;
int_dirs=linspace(min(stims_dir),max(stims_dir),interpfact);
int_grating=interp1(stims_dir,gratings_resp_f2,int_dirs,'spline');
int_plaid=interp1(stims_dir,plaids_resp_f2,int_dirs,'spline');
int_grating=int_grating./max(int_grating(:));
int_plaid=int_plaid./max(int_plaid(:));
plot(int_dirs,int_grating,'color',col1,'linewidth',3)
plot(int_dirs,int_plaid,'color',col1/2,'linewidth',3)
plot([180,180],get(gca,'ylim'),':','color',[0.5,0.5,0.5],'linewidth',2)
xlim([180-90,180+90])
title('Grating and Plaid tuning curves')
ylabel('filter response')
xlabel('direction (?)')
set(gca,'fontsize',12)
title('Grating and Plaid tuning curves')
ylim([0,1.1]);
ylabel('max response (normalized)');
xlabel('direction (?)')
set(gca,'fontsize',12)
% aspect ratio 3 ----------------------------------------------
subplot(3,4,9)
imagesc(Gabor3); colormap('gray'); colorbar;
caxis([-max(abs(Gabor3(:))),+max(abs(Gabor3(:)))]);
axis square;
xlabel('X');
ylabel('Y');
title(['Gabor filter ( aspect ratio = ',num2str(aspect3),' )']);
set(gca,'fontsize',12)
subplot(3,4,10)
imagesc(gratings_resp_dy_f3); colormap(gray); colorbar;
ylabel('timebin #')
xlabel('direction (?)')
title('Grating')
xticks(1:numel(stims_dir)); xtickangle(45);
xticklabels(split(cellstr(num2str(stims_dir))))
axis square
set(gca,'fontsize',12)
subplot(3,4,11)
imagesc(plaids_resp_dy_f3); colormap(gray); colorbar;
ylabel('timebin #')
xlabel('direction (?)')
title('Plaid')
xticks(1:numel(stims_dir)); xtickangle(45);
xticklabels(split(cellstr(num2str(stims_dir))))
axis square
set(gca,'fontsize',12)
subplot(3,4,12)
hold on;
int_dirs=linspace(min(stims_dir),max(stims_dir),interpfact);
int_grating=interp1(stims_dir,gratings_resp_f3,int_dirs,'spline');
int_plaid=interp1(stims_dir,plaids_resp_f3,int_dirs,'spline');
int_grating=int_grating./max(int_grating(:));
int_plaid=int_plaid./max(int_plaid(:));
plot(int_dirs,int_grating,'color',col2,'linewidth',3)
plot(int_dirs,int_plaid,'color',col2/2,'linewidth',3)
plot([180,180],get(gca,'ylim'),':','color',[0.5,0.5,0.5],'linewidth',2)
xlim([180-90,180+90])
title('Grating and Plaid tuning curves')
ylim([0,1.1]);
ylabel('max response (normalized)');
xlabel('direction (?)')
set(gca,'fontsize',12);
suptitle('Blobbyness experiment');
set(gca,'fontsize',12);
% save results
saveas(f0,['blobbyness_visualization'],'jpg')
print(f0,'-depsc','-painters',['blobbyness_visualization','.eps'])