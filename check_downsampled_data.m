%%
function check_downsampled_data(epochT,ops)
t1 = epochT(1);
t2 = epochT(2);
tTotal = t2-t1;
NTbuff = tTotal*ops.fs;
%% load raw data
fid         = fopen(ops.fbinary, 'r'); % open for reading raw data
StartSample = floor(t1*ops.fs);
offset1 = 2*ops.NchanTOT*StartSample; % number of samples to start reading at.
fseek(fid, offset1, 'bof'); % fseek to batch start in raw file
buffa = fread(fid, [ops.NchanTOT NTbuff], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
fclose(fid)
buffa  = buffa(ops.chanMap1,:); % subsample only good channels
%% load LFP data
lfp_file = ops.fproc1;
% load LFP.data binary file
fsize  = get_file_size(lfp_file)/384/2; % size in bytes of raw binary
fidW1  = fopen(lfp_file,   'r'); % open for writing processed data
buffTemp3 = fread(fidW1, [384 fsize],'int16'); % write this batch to binary file
fclose(fidW1);
%%
figure; 
ta = linspace(t1,t2,size(buffa,2));
plot(ta,buffa(100,:),'k');
hold on
buffDown = buffTemp3(100,t1*100+1:floor(t2*100));
% tb = linspace(t1,t2,size(buffDown,2));
tb = t1:1/100:t2;
tb = tb(1:end-1);
plot(tb,buffDown,'r');
end
