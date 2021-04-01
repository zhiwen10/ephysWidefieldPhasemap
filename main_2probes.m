githubDir = 'C:\Users\Steinmetz lab\Documents\git';
%% Add paths
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master')) % path to kilosort folder
addpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\ephysDoubeProbePhasemap')
%% session info
mn = 'ZYE_0020';
td = '2021-03-26';
en = 1;
imec = 1;
imec1 = 0;
probeName = 'p1';
serverRoot = expPath(mn, td, en);
%% 
% save lowpass filtered and dowsampled LFP data to ops.fproc1
ops = configure_file(serverRoot,probeName,imec);
ops.fproc1 = 'E:\ephys\ZYE_0020\2021-03-26\LFP\filtered_lfp_imec1.dat';
ops1 = configure_file(serverRoot,probeName,imec1);
ops1.fproc1 = 'E:\ephys\ZYE_0020\2021-03-26\LFP\filtered_lfp_imec0.dat';
%% check dowsnampled data with raw ephys data and widefield data
epochT = [1200 1250];
%% plot ephys channel phasemap
% change for different probe map profile
[traceEphysMapa,phaseEphysMapa,ampEphysMapa] = ephys_phasemap(ops, epochT);
traceEphysMap1a = reshape(traceEphysMapa,size(traceEphysMapa,1),[]);
phaseEphysMap1a = reshape(phaseEphysMapa,size(phaseEphysMapa,1),[]);
%%
[traceEphysMapb,phaseEphysMapb,ampEphysMapb] = ephys_phasemap(ops1, epochT);
traceEphysMap1b = reshape(traceEphysMapb,size(traceEphysMapb,1),[]);
phaseEphysMap1b = reshape(phaseEphysMapb,size(phaseEphysMapb,1),[]);
%%
traceEphysMap = cat(3,traceEphysMapb,traceEphysMapa);
traceEphysMap1 = reshape(traceEphysMap,size(traceEphysMap,1),[]);
phaseEphysMap = cat(3,phaseEphysMapb,phaseEphysMapa);
phaseEphysMap1 = reshape(phaseEphysMap,size(phaseEphysMap,1),[]);
%%
traceEphysSize = size(traceEphysMap1,1);
%% trace plot
figure; 
hold on
plot(epochT(1):1/100:epochT(1)+(traceEphysSize-1)/100,mean(traceEphysMap1,2)/10+1000,'r')
hold on
plot(epochT(1):1/100:epochT(1)+(traceEphysSize-1)/100,mean(phaseEphysMap1,2)*10+1100,'k')

text(epochT(1)+2,1300,'LFP+WF phase 100Hz');
%%
tEpoch = epochT(1):1/100:epochT(1)+(traceEphysSize-1)/100;
% indx = find(tEpoch>=1121.3 & tEpoch<1121.5);
indx = find(tEpoch>=1015 & tEpoch<1020);
% indx = find(tEpoch>=1124.9 & tEpoch<1125.1);
traceEphysMap2 = traceEphysMap(indx,:,:);
phaseEphysMap2 = phaseEphysMap(indx,:,:);
phaseEphysMap3 = unwrap(reshape(phaseEphysMap2,[],48*16));
phaseEphysMap3 = reshape(phaseEphysMap3,[],48,16);
ampEphysMap2 = ampEphysMap(indx,:,:);
xwf = (1:numel(indx))/numel(indx)*0.8;
%% phasemap video
sizeN = size(phaseEphysMap,1);
immax = max(max(max(phaseEphysMap)));
immin = min(min(min(phaseEphysMap)));
maxTraceEphys = max(max(max(traceEphysMap)));
minTraceEphys = min(min(min(traceEphysMap)));
imH = []; imH2 = [];

ha = figure;
v = VideoWriter(['ephys_trace& phasemap_2probes_ba_' num2str(epochT(1)) '_' num2str(epochT(2)) 's.avi']);
v.FrameRate = 5;
open(v);
for i = 1:sizeN 
    % traceEphysMap
    subplot(1,2,1)
    At = squeeze(traceEphysMap(i,:,:));
    if i == 1
        imHt = imagesc(At);
        axis image; % axis off;
        caxis([minTraceEphys, maxTraceEphys])
        colormap(hsv);
        % set(gca, 'xdir','reverse')
%         frameT = epochT(1)+i/100;
%         h1 = text(1,46,[num2str(frameT) ' s']);
    else
        set(imHt, 'CData', At);
%        frameT = epochT(1)+i/100;
%         h1.String = [num2str(frameT) ' s'];
    end   
    %phaseEphysMap
    subplot(1,2,2)
    A = squeeze(phaseEphysMap(i,:,:));
    if i == 1
        imH = imagesc(A);
        axis image; % axis off;
        caxis([immin, immax])
        colormap(hsv);
        % set(gca, 'xdir','reverse')
        frameT = epochT(1)+i/100;
        h1 = text(1,46,[num2str(frameT) ' s']);
    else
        set(imH, 'CData', A);
        frameT = epochT(1)+i/100;
        h1.String = [num2str(frameT) ' s'];
    end   
    thisFrame = getframe(ha);
    writeVideo(v, thisFrame);
end
close(v);       