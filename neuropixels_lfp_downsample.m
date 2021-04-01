function neuropixels_lfp_downsample(ops)
% weights to combine batches at the edge
w_edge = linspace(0, 1, ops.ntbuff)';
ntb = ops.ntbuff;
NchanTOT = ops.NchanTOT;
NT = ops.NT;
NTbuff = ops.NTbuff;
chanMap1 = ops.chanMap1;
datr_prev = gpuArray.zeros(ntb, ops.Nchan, 'single');
% datr_prev = zeros(ntb, ops.Nchan, 'single');
%%
tDown = [];
buffTemp = [];
Nbatch = ops.Nbatch;
fid         = fopen(ops.fbinary, 'r'); % open for reading raw data
% for ibatch = 1:2
for ibatch = 1:Nbatch
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
    datr = datr(:, chanMap1); % subsample only good channels
    % apply low pass filter at 100Hz
    datr    = gpufilter1(buff, ops, chanMap1); % apply filters and median subtraction
    % datr(ntb + [1:ntb], :) = datr_prev;
    datr(ntb + [1:ntb], :) = w_edge .* datr(ntb + [1:ntb], :) +...
        (1 - w_edge) .* datr_prev;
   
    datr_prev = datr(ntb +NT + [1:ops.ntbuff], :);
    datr    = datr(ntb + (1:NT),:); % remove timepoints used as buffers
    datcpu  = gather(int16(datr')); % convert to int16, and gather on the CPU side
%     buffT = int16(datcpu);
%     count2 = fwrite(fidw, buffT, 'int16');
    %% downsample filtered ephys data to 100Hz
    buffTest1 = double(datcpu);
    batchStart = NT * (ibatch-1)/ops.fs;
    batchEnd = NT*ibatch/ops.fs;
    t1 = NT * (ibatch-1)/ops.fs:1/ops.fs:(NT*ibatch-1)/ops.fs;
    % downsample 300x to 100Hz
    t2 = NT * (ibatch-1)/ops.fs:1/ops.fs*300:(NT*ibatch-1)/ops.fs;
    buffTest2 = interp1(t1,buffTest1',t2);      
    tDown = [tDown t2];
    buffTemp = [buffTemp; buffTest2];
    fprintf([num2str(ibatch) '/' num2str(ops.Nbatch) ' batches \n']);
end
fclose(fid);
% fclose(fidw);
%% interp buffTemp data evenly spaced at 100Hz by using tDown, to get rid of the jitter around the edges of the batches
buffTemp2 = interp1(tDown,buffTemp,tDown(1):1/100:tDown(end));
buffTemp2 = int16(buffTemp2');
% white buffTemp2 to lfp binary file
fidW1        = fopen(ops.fproc1,   'w'); % open for writing processed data
count1 = fwrite(fidW1, buffTemp2, 'int16'); % write this batch to binary file
fclose(fidW1);
end
%% check raw ephys data with filtered & downsampled data between t1 and t2
% t1 = (1-1)*NT/ops.fs;
% t2 = 2*NT/ops.fs;
% tTotal = t2-t1;
% NTbuff = tTotal*ops.fs;
% %%
% % load raw ephys data from t1 to t2
% fid         = fopen(ops.fbinary, 'r'); % open for reading raw data
% StartSample = floor(t1*ops.fs);
% offset1 = 2*NchanTOT*StartSample; % number of samples to start reading at.
% fseek(fid, offset1, 'bof'); % fseek to batch start in raw file
% buffa = fread(fid, [NchanTOT NTbuff], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
% fclose(fid)
% % plot raw ephys data and 
% ta = linspace(t1,t2,size(buffa,2));
% figure
% plot(ta,buffa(100,:),'g');
% hold on
% plot(tDown,buffTemp(:,100),'k')
% hold on
% plot(tDown(1):1/100:tDown(end),buffTemp2(100,:),'r')