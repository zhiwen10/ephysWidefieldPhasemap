function plot_MUA(serverRoot,probeName,imec,ops,epochT)
tt1 = epochT(1); tt2 = epochT(2);
%% load spikes
ksRoot = fullfile(fileparts(getProbeFile(serverRoot, probeName,imec)));
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
%% MUA histogram and plot for 4 shanks
twf = tt1:1/35:tt2;
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
probeTips = chanMap1(30,[1,3,5,7]);
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
    plot(spikeT,spikeHist*3+600,'color',color1{shank})
    hold on; plot(ta,buffa(probeTips(shank),:)/5+500,'color',color1{shank});
    % plot LFP data
    buffDown = buffTemp4(shank,tt1*100+1:floor(tt2*100));
    plot(tb,buffDown+800,'color',color1{shank});
end
text(tt1+2,850,'filtered LFP');
text(tt1+2,700,'MUA');
text(tt1+2,600,'Raw');

end
%% spike_histogram
function [n,x] = spike_histogram(spikeT,binSize)
time_max = max(spikeT);
time_min = min(spikeT);
[n,x] = hist(spikeT, time_min:binSize:time_max);
end