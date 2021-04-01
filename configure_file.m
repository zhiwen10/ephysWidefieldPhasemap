%%
ops.fbinary = getProbeFile(serverRoot, probeName);
% frequency for high pass filtering (150)
ops.fshigh = 50; 
ops.fslow = 100;
ops.NchanTOT    = 385; % total number of channels in your recording
ops.fs = 30000; 
ops.trange = [0 Inf]; % time range to sort
ops.ntbuff = 64;    % samples of symmetrical buffer for whitening and spike detection
%%
ops.NT = 64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 
NchanTOT = ops.NchanTOT; % total number of channels in the raw binary file, including dead, auxiliary etc
bytes       = get_file_size(ops.fbinary); % size in bytes of raw binary
nTimepoints = floor(bytes/NchanTOT/2); % number of total timepoints
%%
ops.tstart  = ceil(ops.trange(1) * ops.fs); % starting timepoint for processing data segment
ops.tend    = min(nTimepoints, ceil(ops.trange(2) * ops.fs)); % ending timepoint
ops.sampsToRead = ops.tend-ops.tstart; % total number of samples to read
ops.twind = ops.tstart * NchanTOT*2; % skip this many bytes at the start
NT = ops.NT;
Nbatch      = ceil(ops.sampsToRead /NT); % number of data batches
ops.Nbatch = Nbatch;
[chanMap1, xc, yc, kcoords, NchanTOTdefault] = loadChanMap(ops.chanMap); % function to load channel map file
ops.chanMap1 = chanMap1;
ops.NchanTOT = getOr(ops, 'NchanTOT', NchanTOTdefault); % if NchanTOT was left empty, then overwrite with the default
ops.Nchan = numel(chanMap1); % total number of good channels that we will spike sort
ops.Nfilt = getOr(ops, 'nfilt_factor', 4) * ops.Nchan; % upper bound on the number of templates we can have
rez.ops         = ops; % memorize ops
ops.NTbuff      = ops.NT + 3*ops.ntbuff; % we need buffers on both sides for filtering
rez.ops.Nbatch = Nbatch;
rez.ops.NTbuff = ops.NTbuff;
rez.ops.chanMap1 = chanMap1;
