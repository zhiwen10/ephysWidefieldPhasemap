githubDir = 'C:\Users\Steinmetz lab\Documents\git';
%% Add paths
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'Pipelines')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master')) % path to kilosort folder
%%
% find the binary file
mn = 'ZYE_0016';
td = '2021-02-25';
en = 1;
imec = 0;
probeName = 'p1';
serverRoot = expPath(mn, td, en);
%%
ops.fbinary = getProbeFile(serverRoot, probeName);
ops.fproc = 'lfp_lpf.dat';
ops.fproc1 = 'lfp_lpf_denoise.dat';
% frequency for high pass filtering (150)
ops.fshigh = 50; 
ops.fslow = 100;
ops.NchanTOT    = 385; % total number of channels in your recording
ops.fs = 30000; 
% ops.ts = 600; %duration of epoch in seconds
ops.trange = [0 Inf]; % time range to sort
ops.ntbuff = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.chanMap = 'C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master\configFiles\NPtype24_hStripe_botRow0_ref1.mat';
%%
tic
ops.nt0 	  = getOr(ops, {'nt0'}, 61); % number of time samples for the templates (has to be <=81 due to GPU shared memory)
ops.nt0min  = getOr(ops, 'nt0min', ceil(20 * ops.nt0/61)); % time sample where the negative peak should be aligned

% NT       = ops.fs*ops.ts; % number of timepoints per batch
NT = 64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 
NchanTOT = ops.NchanTOT; % total number of channels in the raw binary file, including dead, auxiliary etc

bytes       = get_file_size(ops.fbinary); % size in bytes of raw binary
nTimepoints = floor(bytes/NchanTOT/2); % number of total timepoints
ops.tstart  = ceil(ops.trange(1) * ops.fs); % starting timepoint for processing data segment
ops.tend    = min(nTimepoints, ceil(ops.trange(2) * ops.fs)); % ending timepoint
ops.sampsToRead = ops.tend-ops.tstart; % total number of samples to read
ops.twind = ops.tstart * NchanTOT*2; % skip this many bytes at the start

Nbatch      = ceil(ops.sampsToRead /(NT-ops.ntbuff)); % number of data batches
ops.Nbatch = Nbatch;

[chanMap, xc, yc, kcoords, NchanTOTdefault] = loadChanMap(ops.chanMap); % function to load channel map file
ops.NchanTOT = getOr(ops, 'NchanTOT', NchanTOTdefault); % if NchanTOT was left empty, then overwrite with the default

ops.igood = true(size(chanMap));

ops.Nchan = numel(chanMap); % total number of good channels that we will spike sort
ops.Nfilt = getOr(ops, 'nfilt_factor', 4) * ops.Nchan; % upper bound on the number of templates we can have

rez.ops         = ops; % memorize ops

rez.xc = xc; % for historical reasons, make all these copies of the channel coordinates
rez.yc = yc;
rez.xcoords = xc;
rez.ycoords = yc;
% rez.connected   = connected;
rez.ops.chanMap = chanMap;
rez.ops.kcoords = kcoords;


NTbuff      = NT + 3*ops.ntbuff; % we need buffers on both sides for filtering

rez.ops.Nbatch = Nbatch;
rez.ops.NTbuff = NTbuff;
rez.ops.chanMap = chanMap;
% weights to combine batches at the edge
w_edge = linspace(0, 1, ops.ntbuff)';
ntb = ops.ntbuff;
% datr_prev = gpuArray.zeros(ntb, ops.Nchan, 'single');
datr_prev = zeros(ntb, ops.Nchan, 'single');
fid         = fopen(ops.fbinary, 'r'); % open for reading raw data
fidW        = fopen(ops.fproc,   'w'); % open for writing processed data
fidW1        = fopen(ops.fproc1,   'w'); % open for writing processed data
%%
syncTL = loadAlign(serverRoot, 'tl');
syncProbe = loadAlign(serverRoot, probeName);
if size(syncProbe,1)-size(syncTL,1) ~= 0
    syncTL = syncTL(1:size(syncProbe,1),1);
end
wideExp = readNPY(fullfile(serverRoot,'widefieldExposure.raw.npy'));
ts = readNPY(fullfile(serverRoot,'widefieldExposure.timestamps_Timeline.npy'));
t = tsToT(ts,size(wideExp,1));
[~, frTimes] = schmittTimes(t, wideExp, [1 2]); 
frTimes1 = interp1(syncTL, syncProbe, frTimes);  
%%
ibatch1 = 501;
ibatch2 = 510;
%%
tDown = [];
for ibatch = ibatch1:ibatch2
    % we'll create a binary file of batches of NT samples, which overlap consecutively on ops.ntbuff samples
    % in addition to that, we'll read another ops.ntbuff samples from before and after, to have as buffers for filtering
    offset = max(0, ops.twind + 2*NchanTOT*(NT * (ibatch-1) - ntb)); % number of samples to start reading at.
    fseek(fid, offset, 'bof'); % fseek to batch start in raw file

    buff = fread(fid, [NchanTOT NTbuff], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
    if isempty(buff)
        break; % this shouldn't really happen, unless we counted data batches wrong
    end
    nsampcurr = size(buff,2); % how many time samples the current batch has
    if nsampcurr<NTbuff
        buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr); % pad with zeros, if this is the last batch
    end
    if offset==0
        bpad = repmat(buff(:,1), 1, ntb);
        buff = cat(2, bpad, buff(:, 1:NTbuff-ntb)); % The very first batch has no pre-buffer, and has to be treated separately
    end
    datr = gpuArray(buff);
    datr = single(datr');
    datr = datr(:, chanMap); % subsample only good channels
    datr    = gpufilter1(buff, ops, chanMap); % apply filters and median subtraction
    datr(ntb + [1:ntb], :) = datr_prev;
    datr(ntb + [1:ntb], :) = w_edge .* datr(ntb + [1:ntb], :) +...
        (1 - w_edge) .* datr_prev;
   
    datr_prev = datr(ntb +NT + [1:ops.ntbuff], :);
    datr    = datr(ntb + (1:NT),:); % remove timepoints used as buffers
    datcpu  = gather(int16(datr')); % convert to int16, and gather on the CPU side
    count = fwrite(fidW, datcpu, 'int16'); % write this batch to binary file
    if count~=numel(datcpu)
        error('Error writing batch %g to %s. Check available disk space.',ibatch,ops.fproc);
    end
    %%
    % replace 80 ephys samples near wfExposure (LED artifacts) with values
    % interpreted from two joint points
    batchStart = NT * (ibatch-1)/ops.fs;
    batchEnd = NT*ibatch/ops.fs;
    t1 = NT * (ibatch-1)/ops.fs:1/ops.fs:(NT*ibatch-1)/ops.fs;
    frTimes2 = frTimes1(frTimes1>batchStart & frTimes1<batchEnd)-0.002;
    buffTest1 = double(datcpu);
    for i = 1:numel(frTimes2)
        [~,tTemp1] = min(abs(t1-frTimes2(i)));    
        if tTemp1+80 <= size(t1,2)
            tt1 = t1(tTemp1); tTrace1 = buffTest1(:,tTemp1);
            tt2 = t1(tTemp1+80); tTrace2 = buffTest1(:,tTemp1+80);
            tt = [tt1,tt2]; tTrace = [tTrace1, tTrace2];
            buffTemp = interp1(tt,tTrace',tt1:1/ops.fs:tt2); 
            buffTest1(:,tTemp1:tTemp1+80) = buffTemp';
        end
    end
    t2 = NT * (ibatch-1)/ops.fs:1/ops.fs*300:(NT*ibatch-1)/ops.fs;
    buffTest2 = interp1(t1,buffTest1',t2);    
    buffTest2 = int16(buffTest2');    
    count1 = fwrite(fidW1, buffTest2, 'int16'); % write this batch to binary file
    tDown = [tDown t2];
end
%%
fclose(fidW); % close the files
fclose(fid);
fclose(fidW1);
rez.temp.Nbatch = Nbatch;
toc
%%
batchStart = NT * (ibatch1-1)/ops.fs;
batchEnd = (NT*ibatch2-1)/ops.fs;
t1 = NT * (ibatch1-1)/ops.fs:1/ops.fs:(NT*ibatch2-1)/ops.fs;
t2 = downsample(t1,299);
%%
samples = get_file_size('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\np24lfp\lfp.dat')/(384*2); % size in bytes of raw binary
fid2 = fopen('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\np24lfp\lfp.dat', 'r'); % open for reading raw data
buff1 = fread(fid2, [384 samples], '*int16'); 
fclose(fid2);
%%
samples1 = get_file_size('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\np24lfp\lfp_lpf.dat')/(384*2); % size in bytes of raw binary
fid3 = fopen('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\np24lfp\lfp_lpf.dat', 'r'); % open for reading raw data
buff2 = fread(fid3, [384 samples1], '*int16'); 
fclose(fid3);
%%
samples2 = get_file_size('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\np24lfp\lfp_lpf_denoise.dat')/(384*2); % size in bytes of raw binary
fid4 = fopen('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\np24lfp\lfp_lpf_denoise.dat', 'r'); % open for reading raw data
buff3 = fread(fid4, [384 samples2], '*int16'); 
fclose(fid4);
%%
buffTest = double(buff1(100,:));
%%
frTimes2 = frTimes1(frTimes1>batchStart & frTimes1<batchEnd)-0.002;
scatter(frTimes2,2000*ones(size(frTimes2,1),1),'b')
%%
buffTest1 = buffTest;
for i = 1:numel(frTimes2)
    [~,tTemp1] = min(abs(t1-frTimes2(i)));    
    tt1 = t1(tTemp1); tTrace1 = buffTest(tTemp1);
    tt2 = t1(tTemp1+80); tTrace2 = buffTest(tTemp1+80);
    buffTest1(tTemp1:tTemp1+80) = interp1([tt1,tt2],[tTrace1,tTrace2],tt1:1/ops.fs:tt2); 
end
figure;plot(t1,buffTest,'k')
hold on
plot(t1,buffTest1+100,'r')
%%
figure;plot(t1,buff1(100,:),'k')
hold on
plot(t1,buff2(100,:)+50,'r')
hold on
plot(tDown,buff3(100,:)+100,'g')
legend('raw','lpf','lpf_denoise')