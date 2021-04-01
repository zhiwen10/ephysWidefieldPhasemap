function plot_MUA_wfpixel(serverRoot,probeName,ops,epochT,point)
tt1 = epochT(1); tt2 = epochT(2);
%% load spikes
ksRoot = fullfile(fileparts(getProbeFile(serverRoot, probeName)));
sp = loadKSdir(ksRoot);
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(ksRoot);
%% shank sites
if contains(ops.chanMap,'NPtype24_hStripe')
    shanks(:,1) = find(sp.xcoords<-200);
    shanks(:,2) = find(sp.xcoords<100 & sp.xcoords>-200);
    shanks(:,3) = find(sp.xcoords<300 & sp.xcoords>100);
    shanks(:,4) = find(sp.xcoords>300);
elseif contains(ops.chanMap,'NPtype24_doubleLengthStripe')
    shanks(:,1) = find(sp.xcoords<100);
    shanks(:,2) = find(sp.xcoords<300 & sp.xcoords>200);
    shanks(:,3) = find(sp.xcoords<600 & sp.xcoords>400);
    shanks(:,4) = find(sp.xcoords>700);
end
%%
%plot PCA components and traces for blue and purple channels
corrPath = fullfile(serverRoot, 'corr', 'svdTemporalComponents_corr.npy');
if ~exist(corrPath, 'file')
    colors = {'blue', 'violet'};
    computeWidefieldTimestamps(serverRoot, colors); % preprocess video
    nSV = 200;
    [U, V, t, mimg] = hemoCorrect(serverRoot, nSV); % process hemodynamic correction
else
    nSV = 200;
    [U, V, t, mimg] = loadUVt(serverRoot, nSV);
end
if length(t) > size(V,2)
  t = t(1:size(V,2));
end
%%
% pixelCorrelationViewerSVD(U,V)
dV = [zeros(size(V,1),1) diff(V,[],2)];
ddV = [zeros(size(dV,1),1) diff(dV,[],2)];
%%
tlFile = fullfile(serverRoot, 'blue\meanImage.npy');
meanImage = readNPY(tlFile);
%%
% point(4,:) = [235,274]; point(3,:) = [225,318]; point(2,:) = [217,370]; point(1,:) = [209,397];
% for kk = 1:4
%     hold on
%     scatter(point(kk,2),point(kk,1),12,'k','LineWidth',2);
% end
%% realign the time of wf trace to ephys activity
syncTL = loadAlign(serverRoot, 'tl');
% syncProbe = loadAlign(serverRoot, [probeName '_imec' num2str(imec)]);
syncProbe = loadAlign(serverRoot, probeName);
if size(syncProbe,1)-size(syncTL,1) ~= 0
    % syncProbe = [syncProbe; zeros(size(syncTL,1)-size(syncProbe,1), 1)];
    syncTL = syncTL(1:size(syncProbe,1),1);
end
tAll1 = interp1(syncTL, syncProbe, t);  
%% filter wf trace at 2-8Hz with filtfilt
for i  = 1:4
    px_mean = meanImage(point(i,1),point(i,2));
    px1 = squeeze(U(point(i,1),point(i,2),:))'*dV;
    px(i,:) = px1/px_mean;
end
Fs = 35;
[f1,f2] = butter(2, [2 8]/(Fs/2), 'bandpass');
pxTemp = filtfilt(f1,f2,px');
px = pxTemp';
%% MUA histogram and plot for 4 shanks
tAll2 = tAll1(not(isnan(tAll1)));
px2 = px(:,not(isnan(tAll1)));
twf = tt1:1/35:tt2;
pxTemp = interp1(tAll2,px2',twf);
pxTemp = pxTemp';
binSize = 0.025;
% figure
color1 = {'k','r','g','c'};
for shank = 1:4
    clear spikeT spikeHist 
    incl1 = (spikeAmps>20 & ismember(spikeSites,...
            shanks(:,shank)));%& spikeTimes>(3000/35) & spikeTimes<(3500/35));
    s1.spikeTimes = spikeTimes(incl1);
    s1.spikeDepths = spikeDepths(incl1);
    s2{shank} = s1;
    [spikeHist,spikeT] = spike_histogram(s2{shank}.spikeTimes,binSize);
%     plot(spikeTimes(incl1), spikeDepths(incl1),'.','color',color1{shank})
%     hold on
    sIndex  = find(spikeT>=tt1 & spikeT <=tt2);
    spikeT = spikeT(sIndex);
    spikeHist = spikeHist(sIndex);
%     plot(spikeT,spikeHist+1000,'color',color1{shank})
%     hold on
%     plot(tAll1, px(shank,:)*5000+800,'color',color1{shank}); 
%     hold on
%     plot(t2,lfp(probeTip(shank),:)/20+1200,'color',color1{shank});
end
%% load raw ephys data
NchanTOT = ops.NchanTOT;
fs = ops.fs;
tTotal = tt2-tt1;
NTbuff = tTotal*fs;
fid         = fopen(ops.fbinary, 'r'); % open for reading raw data
StartSample = floor(tt1*fs);
offset1 = 2*NchanTOT*StartSample; % number of samples to start reading at.
fseek(fid, offset1, 'bof'); % fseek to batch start in raw file
buffa = fread(fid, [NchanTOT NTbuff], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
fclose(fid)
Map = ops.chanMap;
% [chanMap, xc, yc, kcoords, NchanTOTdefault] = loadChanMap(Map); % function to load channel map file
buffa  = buffa(ops.chanMap1,:); % subsample only good channels
chanMap1 = chanMapReorder(Map);
probeTips = chanMap1(20,[1,3,5,7]);
% probeTips = chanMap1(20,:);
ta = linspace(tt1,tt2,size(buffa,2));
%% load filtered LFP data and filter around 2-8Hz
fsize  = get_file_size(ops.fproc1)/384/2; % size in bytes of raw binary
fidW1  = fopen(ops.fproc1,   'r'); % open for writing processed data
buffTemp3 = fread(fidW1, [384 fsize],'int16'); % write this batch to binary file
fclose(fidW1);
[f1,f2] = butter(3, [2 8]/(100/2), 'bandpass');
buffTemp4 = filtfilt(f1,f2,buffTemp3');
buffTemp4 = buffTemp4';
tb = tt1:1/100:tt2;
tb = tb(1:end-1);
%%
figure;
for shank = 1:4
    clear spikeT spikeHist 
    % plot MUA activity
    incl1 = (spikeAmps>20 & ismember(spikeSites,...
            shanks(:,shank)));%& spikeTimes>(3000/35) & spikeTimes<(3500/35));
    s1.spikeTimes = spikeTimes(incl1);
    s1.spikeDepths = spikeDepths(incl1);
    s2{shank} = s1;
    [spikeHist,spikeT] = spike_histogram(s2{shank}.spikeTimes,binSize);
%     plot(spikeTimes(incl1), spikeDepths(incl1),'.','color',color1{shank})
%     hold on
    sIndex  = find(spikeT>=tt1 & spikeT <=tt2);
    spikeT = spikeT(sIndex);
    spikeHist = spikeHist(sIndex);
    plot(spikeT,spikeHist+600,'color',color1{shank})
    hold on
    % plot raw ephys trace
    plot(twf, pxTemp(shank,:)*10000+200,'color',color1{shank}); 
    % plot plot raw data at the 4 shank sites
    hold on; plot(ta,buffa(probeTips(shank),:)/20+500,'color',color1{shank});
    % plot LFP data
    buffDown = buffTemp4(shank,tt1*100+1:floor(tt2*100));
    plot(tb,buffDown+800,'color',color1{shank});
end
text(tt1+2,850,'filtered LFP');
text(tt1+2,700,'MUA');
text(tt1+2,500,'Raw');
text(tt1+2,300,'WideField');
end
%% spike_histogram
function [n,x] = spike_histogram(spikeT,binSize)
time_max = max(spikeT);
time_min = min(spikeT);
[n,x] = hist(spikeT, time_min:binSize:time_max);
end