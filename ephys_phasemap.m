function [traceEphysMap,phaseEphysMap,ampEphysMap] = ephys_phasemap(ops, epochT)
%% load downsampled LFP data
lfp_file = ops.fproc1;
fsize  = get_file_size(lfp_file)/384/2; % size in bytes of raw binary
fidW1  = fopen(lfp_file, 'r'); % open for writing processed data
lfp = fread(fidW1, [384 fsize],'int16'); % write this batch to binary file
fclose(fidW1);
%%
chanMap1 = chanMapReorder(ops.chanMap);
%%
for j = 1:384
    [row(j),col(j)] = ind2sub(size(chanMap1),find(chanMap1 == j));
    lfpmat(row(j),col(j),:) = lfp(j,:);
end
% for i = 1:48
%     for j = 1:8
%         mapIndex = chanMap1(i,j);
%         lfpmat(i,j,:) = lfp(mapIndex,:);
%     end
% end   
%%
meanTrace = reshape(lfpmat, 384, size(lfpmat,3));
%% 
Fs = 100;
t2 = 0:1/Fs:size(lfpmat,3)/Fs;
t2 = t2(1:end-1);
%%
[min1,indx1] = min(abs(t2-epochT(1)));
[min2,indx2] = min(abs(t2-epochT(2)));
% [min1,indx1] = min(abs(t2-898.5));
% [min2,indx2] = min(abs(t2-899.1));
%% padding 100 sample before and after epoch for interp1 at tTest1 and tTest2
EphysIndx = indx1:indx2;
EphysIndx2 = indx1-100:indx2+100;
meanTrace = meanTrace(:,EphysIndx2);
%% interp1 at exact time of tTest1 and tTest2
rate = 1/Fs;
tq1 = epochT(1):rate:epochT(2);
meanTrace1 = interp1(t2(EphysIndx2),meanTrace',tq1);
meanTrace1 = meanTrace1';
%% filter 2-8Hz
traceEphys = meanTrace1-repmat(mean(meanTrace1,2),1,size(meanTrace1,2));
% filter and hilbert transform work on each column
traceEphys = traceEphys';
[f1,f2] = butter(2, [2 8]/(Fs/2), 'bandpass');
traceEphys = filtfilt(f1,f2,traceEphys);
traceEphysMap = reshape(traceEphys,size(traceEphys,1),48,8); 
% traceEphysMap = reshape(traceEphys,size(traceEphys,1),96,4); 
traceHilbert =hilbert(traceEphys);
tracePhase1 = angle(traceHilbert);
traceAmp1 = abs(traceHilbert);
phaseEphysMap = reshape(tracePhase1,size(tracePhase1,1),48,8); 
ampEphysMap = reshape(traceAmp1,size(traceAmp1,1),48,8); 
% phaseEphysMap = reshape(tracePhase1,size(tracePhase1,1),96,4); 
end
%% ephys_phasemap video
% figure; 
% for i = 1:size(tracePhaseEphys,1)
% imagesc(squeeze(tracePhaseEphys(i,:,:)));
% axis image
% colormap(hsv); 
% frameT = i/350;
% text(1,46,[num2str(frameT) ' s']);
% % axis image; axis off;
% pause(0.2)
% % set(gca, 'ydir','reverse')
% end
% 
% %%
% tracePhase1 = tracePhaseEphys(indx1:indx2,:,:);
% sizeN = size(tracePhase1,1);
% figure;
% immax = max(max(max(tracePhase1)));
% immin = min(min(min(tracePhase1)));
% imH = [];
% v = VideoWriter(['ephys_phasemap_' num2str(indx1) '_' num2str(indx2) '.avi']);
% v.FrameRate = 5;
% open(v);
% for i = 1:sizeN 
%     A = squeeze(tracePhase1(i,:,:));
%     if i == 1
%         imH = imagesc(A);
%         axis image; axis off;
%         caxis([immin, immax])
%         colormap(hsv);
%         
%         frameT = (indx1+i)/100;
%         h1 = text(1,46,[num2str(frameT) ' s']);
%     else
%         set(imH, 'CData', A);
%         frameT = (indx1+i)/100;
%         h1.String = [num2str(frameT) ' s'];
%     end 
%     thisFrame = getframe(gca);
%     writeVideo(v, thisFrame);
% end
% close(v);       