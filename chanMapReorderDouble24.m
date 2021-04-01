% chanMap = 'C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master\configFiles\NPtype24_hStripe_botRow0_ref1.mat';
function chanMap1 = chanMapReorderDouble24(Map)
%% reorder probes to 4 shanks
[chanMap, xc, yc, kcoords, NchanTOTdefault] = loadChanMap(Map); % function to load channel map file
%%
xc1 = xc;
yc1 = floor((yc)/15)+1;
mapReal = [xc1,yc1];
%% find chanMap order to reorder ephys data to 2d map
map1(1:48,1) = zeros(48,1); map1(1:48,2) = 1:48;
map2 = [];
for i = 1:8
    mapTemp = map1;
    if mod(i,2) ==0
        mapTemp(:,1) = mapTemp(:,1)+32+250*(i/2-1);
    else
        mapTemp(:,1) = mapTemp(:,1)+250*(i-1)/2;
    end
    map2 = [map2;mapTemp];
end
[a,b] = ismember(map2,mapReal,'row');
chanMap1 = reshape(b,[48,8]);