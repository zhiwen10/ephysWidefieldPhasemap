githubDir = 'C:\Users\Steinmetz lab\Documents\git';
%% Add paths
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'Pipelines')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master')) % path to kilosort folder
addpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\ephysWidefieldPhasemap')
addpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\pinwheelDetection')
%% session info
mn = 'ZYE_0016';
td = '2021-02-25';
en = 1;
imec = 0;
probeName = 'p1';
serverRoot = expPath(mn, td, en);
%% 
% save lowpass filtered and dowsampled LFP data to ops.fproc1
ops.fproc1 = 'E:\ephys\ZYE_0016\2021-02-25\lfp_filtered_downSample_all.dat';
ops.chanMap = 'C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master\configFiles\NPtype24_hStripe_botRow0_ref1.mat';
% ops.chanMap = 'C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master\configFiles\NPtype24_doubleLengthStripe_botRow0_ref0.mat';
run(fullfile('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\ephysWidefieldPhasemap', 'configure_file'))
%% generate filtered and downsampled LFP data
% neuropixels_lfp_downsample(ops);
%% check downsampled data with raw ephys data
epochT = [1100 1150];
% epochT = [1117 1125];
check_downsampled_data(epochT,ops);
%% check dowsnampled data with raw ephys data and widefield data
% tt1 = 1000; tt2 = 1020;
point(4,:) = [217,312]; point(3,:) = [218,326]; point(2,:) = [215,340]; point(1,:) = [214,353];
plot_MUA_wfpixel(serverRoot,probeName,ops,epochT,point);
%% time points
% tTest = [514.6, 515.3];
% tTest = [898.5 899.1];
syncTL = loadAlign(serverRoot, 'tl');
syncProbe = loadAlign(serverRoot, probeName);
epochT1 = interp1(syncProbe, syncTL, epochT); 
%% plot ephys channel phasemap
% change for different probe map profile
[traceEphysMap,phaseEphysMap,ampEphysMap] = ephys_phasemap(ops, epochT);
traceEphysMap1 = reshape(traceEphysMap,size(traceEphysMap,1),[]);
phaseEphysMap1 = reshape(phaseEphysMap,size(phaseEphysMap,1),[]);
%% plot widefiled phasemap
[traceWF,phaseWF] = get_widefield_phasemap(serverRoot,probeName,epochT1);
%%
traceEphysSize = size(traceEphysMap1,1);
traceWFSize = size(traceWF,1);
%% trace plot
% figure; 
hold on
plot(epochT(1):1/100:epochT(1)+(traceEphysSize-1)/100,mean(traceEphysMap1,2)/10+1000,'r')
hold on

plot(epochT(1):1/100:epochT(1)+(traceWFSize-1)/100,squeeze(traceWF(:,point(1,1),point(1,2)))/10+1000,'k')

hold on
plot(epochT(1):1/100:epochT(1)+(traceEphysSize-1)/100,mean(phaseEphysMap1,2)*10+1200,'r')
hold on
plot(epochT(1):1/100:epochT(1)+(traceWFSize-1)/100,squeeze(phaseWF(:,point(1,1),point(1,2)))*10+1200,'k')

text(epochT(1)+2,1300,'LFP+WF phase 100Hz');
%%
tEpoch = epochT(1):1/100:epochT(1)+(traceEphysSize-1)/100;
% indx = find(tEpoch>=1121.3 & tEpoch<1121.5);
indx = find(tEpoch>=1120.9 & tEpoch<1121.9);
% indx = find(tEpoch>=1124.9 & tEpoch<1125.1);
traceEphysMap2 = traceEphysMap(indx,:,:);
traceEphysMap3 = traceEphysMap2-traceEphysMap2(:,30,8);
phaseEphysMap2 = phaseEphysMap(indx,:,:);
phaseEphysMap3 = unwrap(reshape(phaseEphysMap2,[],48*8));
phaseEphysMap3 = reshape(phaseEphysMap3,[],48,8);
ampEphysMap2 = ampEphysMap(indx,:,:);
traceWFMap2 = traceWF(indx,:,:);
phaseWFMap2 = phaseWF(indx,:,:);
xwf = (1:numel(indx))/numel(indx)*0.8;

%%
figure
subplot(1,3,1)
for i = 1:48
    for j = 1:8
    plot(j+xwf-0.5, i+0.5*(traceEphysMap2(:,i,j)-traceEphysMap2(1,i,j))*0.002,'k');
    hold on;
    xline(j+xwf(50)-0.5,'r');
    hold on;
    end
end
xlim([0 9])
ylim([0 50])
set(gca, 'ydir','reverse')
title('filtered_LFP')
subplot(1,3,2)
for i = 1:48
    for j = 1:8
    plot(j+xwf-0.5, i+0.5*phaseEphysMap3(:,i,j)*0.2,'k');
    hold on;
    xline(j+xwf(50)-0.5,'r');
    hold on;
    end
end
xlim([0 9])
ylim([0 50])
set(gca, 'ydir','reverse')
title('phaseMap')
subplot(1,3,3)
for i = 1:48
    for j = 1:8
        y = ampEphysMap2(:,i,j)-ampEphysMap2(1,i,j);
        plot(j+xwf-0.5, i+0.5*y*0.003,'k');
        hold on;
        xline(j+xwf(50)-0.5,'r');
        hold on;
    end
end
xlim([0 9])
ylim([0 50])
set(gca, 'ydir','reverse')
title('ampMap')
%%
figure;
frameN = size(phaseEphysMap3,1);
for i = 1:frameN
imagesc(squeeze(phaseEphysMap3(i,:,:)));
% set(gca, 'xdir','reverse')
axis image
colormap(hsv); 
frameT = i/100;
text(1,46,[num2str(frameT) ' s']);
axis image; % axis off;
saveas(gcf,['frame' num2str(i) '.png'])
end
%%
%%
figure;
frameN = size(phaseEphysMap3,1);
for i = 1:frameN
imagesc(squeeze(phaseEphysMap3(i,:,:)));
% set(gca, 'xdir','reverse')
axis image
colormap(hsv); 
frameT = i/100;
text(1,46,[num2str(frameT) ' s']);
axis image; % axis off;
saveas(gcf,['frame' num2str(i) '.png'])
end
%% phasemap video
sizeN = min(size(phaseEphysMap,1),size(phaseWF,1));

immax = max(max(max(phaseEphysMap)));
immin = min(min(min(phaseEphysMap)));
immax2 = max(max(max(phaseWF)));
immin2 = min(min(min(phaseWF)));
maxTraceEphys = max(max(max(traceEphysMap)));
minTraceEphys = min(min(min(traceEphysMap)));
maxTraceWF = max(max(max(traceWF)));
minTraceWF = min(min(min(traceWF)));
imH = []; imH2 = [];

ha = figure;
v = VideoWriter(['ephys_trace& phasemap_' num2str(epochT(1)) '_' num2str(epochT(2)) 's.avi']);
v.FrameRate = 5;
open(v);
for i = 1:sizeN 
    % traceEphysMap
    subplot(2,2,1)
    At = squeeze(traceEphysMap(i,:,:));
    if i == 1
        imHt = imagesc(At);
        axis image; % axis off;
        caxis([minTraceEphys, maxTraceEphys])
        colormap(hsv);
        % set(gca, 'xdir','reverse')
        % frameT = epochT(1)+i/100;
        % h1 = text(1,46,[num2str(frameT) ' s']);
    else
        set(imHt, 'CData', At);
       % frameT = epochT(1)+i/100;
       %  h1.String = [num2str(frameT) ' s'];
    end   
    %phaseEphysMap
    subplot(2,2,2)
    A = squeeze(phaseEphysMap(i,:,:));
    if i == 1
        imH = imagesc(A);
        axis image; % axis off;
        caxis([immin, immax])
        colormap(hsv);
        % set(gca, 'xdir','reverse')
        frameT = epochT(1)+i/100;
        % h1 = text(1,46,[num2str(frameT) ' s']);
    else
        set(imH, 'CData', A);
        frameT = epochT(1)+i/100;
       %  h1.String = [num2str(frameT) ' s'];
    end
    % traceWF
    subplot(2,2,3)
    Bt = squeeze(traceWF(i,:,:));
    if i == 1
        imH2t = imagesc(Bt);
        axis image; % axis off;
        caxis([minTraceWF, maxTraceWF])
        colormap(hsv);
        
        % frameT = epochT(1)+i/100;
        % h2 = text(1,46,[num2str(frameT) ' s']);
    else
        set(imH2t, 'CData', Bt);
        % frameT = epochT(1)+i/100;
        % h2.String = [num2str(frameT) ' s'];
    end 
    
    subplot(2,2,4)
    B = squeeze(phaseWF(i,:,:));
    if i == 1
        imH2 = imagesc(B);
        axis image; % axis off;
        caxis([immin2, immax2])
        colormap(hsv);
        
        frameT = epochT(1)+i/100;
        h2 = text(1,46,[num2str(frameT) ' s']);
    else
        set(imH2, 'CData', B);
        frameT = epochT(1)+i/100;
        h2.String = [num2str(frameT) ' s'];
    end 
    
    thisFrame = getframe(ha);
    writeVideo(v, thisFrame);
end
close(v);       
%%
figure
for i = 1:48
    for j = 1:8
        plot(j+xwf-0.5, i+0.5*(traceEphysMap2(:,i,j))*0.01,'k');
        % plot(j+xwf-0.5, i+0.5*(traceEphysMap2(:,i,j)-traceEphysMap2(:,30,8))*0.01,'k');
        hold on;
        xline(j+xwf(10)-0.5,'r');
        hold on;
    end
end
xlim([0 9])
ylim([0 50])
set(gca, 'ydir','reverse')
title('filtered_LFP')
%%
%% phasemap video2
sizeN = min(size(traceEphysMap3,1),size(phaseWFMap2,1));

immax = max(max(max(phaseEphysMap2)));
immin = min(min(min(phaseEphysMap2)));
immax2 = max(max(max(phaseWFMap2)));
immin2 = min(min(min(phaseWFMap2)));
maxTraceEphys = max(max(max(traceEphysMap3)));
minTraceEphys = min(min(min(traceEphysMap3)));
maxTraceWF = max(max(max(traceWFMap2)));
minTraceWF = min(min(min(traceWFMap2)));
imH = []; imH2 = [];
epochT1 = 1120.9;

ha = figure;
v = VideoWriter(['ephys_trace& phasemap_1120.9-1121.9s-2.avi']);
v.FrameRate = 5;
open(v);
for i = 1:sizeN 
    % traceEphysMap
    subplot(2,2,1)
    At = squeeze(traceEphysMap3(i,:,:));
    if i == 1
        imHt = imagesc(At);
        axis image; % axis off;
        caxis([minTraceEphys, maxTraceEphys])
        colormap(hsv);
        % set(gca, 'xdir','reverse')
        % frameT = epochT(1)+i/100;
        % h1 = text(1,46,[num2str(frameT) ' s']);
    else
        set(imHt, 'CData', At);
       % frameT = epochT(1)+i/100;
       %  h1.String = [num2str(frameT) ' s'];
    end   
    %phaseEphysMap
    subplot(2,2,2)
    A = squeeze(phaseEphysMap2(i,:,:));
    if i == 1
        imH = imagesc(A);
        axis image; % axis off;
        caxis([immin, immax])
        colormap(hsv);
        % set(gca, 'xdir','reverse')
        frameT = epochT1(1)+i/100;
        % h1 = text(1,46,[num2str(frameT) ' s']);
    else
        set(imH, 'CData', A);
        frameT = epochT1(1)+i/100;
       %  h1.String = [num2str(frameT) ' s'];
    end
    % traceWF
    subplot(2,2,3)
    Bt = squeeze(traceWFMap2(i,:,:));
    if i == 1
        imH2t = imagesc(Bt);
        axis image; % axis off;
        caxis([minTraceWF, maxTraceWF])
        colormap(hsv);
        
        % frameT = epochT(1)+i/100;
        % h2 = text(1,46,[num2str(frameT) ' s']);
    else
        set(imH2t, 'CData', Bt);
        % frameT = epochT(1)+i/100;
        % h2.String = [num2str(frameT) ' s'];
    end 
    
    subplot(2,2,4)
    B = squeeze(phaseWFMap2(i,:,:));
    if i == 1
        imH2 = imagesc(B);
        axis image; % axis off;
        caxis([immin2, immax2])
        colormap(hsv);
        
        frameT = epochT1(1)+i/100;
        h2 = text(1,46,[num2str(frameT) ' s']);
    else
        set(imH2, 'CData', B);
        frameT = epochT1(1)+i/100;
        h2.String = [num2str(frameT) ' s'];
    end 
    
    thisFrame = getframe(ha);
    writeVideo(v, thisFrame);
end
close(v);       