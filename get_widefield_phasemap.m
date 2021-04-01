
function [traceWF,phaseWF] = get_widefield_phasemap(serverRoot,probeName,epochT1)
nSV = 200;
[U, V, t, mimg] = loadUVt(serverRoot, nSV);
if length(t) > size(V,2)
  t = t(1:size(V,2));
end
dV = [zeros(size(V,1),1) diff(V,[],2)];
ddV = [zeros(size(dV,1),1) diff(dV,[],2)];
%%
syncTL = loadAlign(serverRoot, 'tl');
syncProbe = loadAlign(serverRoot, probeName);
tAll1 = t;
% tAll1 = interp1(syncTL, syncProbe, t);  
%%
[min1,indx1] = min(abs(tAll1-epochT1(1)));
[min2,indx2] = min(abs(tAll1-epochT1(2)));
winSamps2 = [indx1:indx2];
tSamps2 = tAll1(winSamps2);
winSamps3 = [indx1-35:indx2+35];
tSamps3 = tAll1(winSamps3);
%%
periEventV = dV(:,winSamps3);
Ur = reshape(U, size(U,1)*size(U,2), size(U,3));
%% trace re-construction
meanTrace = Ur*periEventV;
Fs = 100;
rate = 1/Fs;
tq = epochT1(1):rate:epochT1(2);
meanTrace1 = interp1(tSamps3,meanTrace',tq);
meanTrace1 = meanTrace1';
%% filter 2-8Hz
traceTemp = meanTrace1-repmat(mean(meanTrace1,2),1,size(meanTrace1,2));
% filter and hilbert transform work on each column
traceTemp = traceTemp';
[f1,f2] = butter(2, [2 8]/(Fs/2), 'bandpass');
traceTemp = filter(f1,f2,traceTemp);
traceWF = reshape(traceTemp,size(traceTemp,1),size(U,1), size(U,2));
traceHilbert =hilbert(traceTemp);
tracePhase = angle(traceHilbert);
phaseWF = reshape(tracePhase,size(tracePhase,1),size(U,1), size(U,2));
end
%%
% figure; 
% for i = 1:size(tracePhase,1)
% imagesc(squeeze(tracePhase(i,:,:)));
% axis image
% colormap(hsv); 
% frameT = i/350;
% text(1,46,[num2str(frameT) ' s']);
% % axis image; axis off;
% pause(0.05)
% % set(gca, 'ydir','reverse')
% end
