%run.m - example of the run script for the laser speckle contrast analysis
% Author: DD Postnov, PhD
% BOAS lab, Boston University
% email address: dpostnov@bu.edu
% Last revision: 3-March-2018

%% Read and process at the same time - the most standard procedure
clear all
close all
folderName='D:\Evren\trial 1 left'; % folder with raw files
burstStart=1; % first burst to read
burstN=[]; % number of bursts to read, [] for all bursts
kernelSize=7; % kernel size for LSCI processing. 7 or 5 for sLSCI, 25 for tLSCI
lscType='sLSCI'; % type of contrast analysis
procTypeLSCI='gpu'; %processing type for contrast analysis
dsType='none'; %data downsampling (none or kernel)
procType='cluster'; % single thread or parallel processing of multiple raw files

[LSCI,time,imgsPerBurst]=getLSCIFolderIMG(procType,folderName,burstStart,burstN,lscType,kernelSize,procTypeLSCI,dsType);
% substract first time value (in case if no loop over burst files is used)
time=time-time(1);
% save results
save([folderName,'\LSCI','.mat'],'LSCI','time','lscType','dsType','imgsPerBurst','-v7.3');

%% Convert contrast to flow

% assuming that flow is 1/K^2
flowIndex=1./(LSCI.*LSCI);
figure
plot(time,squeeze(mean(mean(LSCI,1),2)));
xlabel('Time,s')
ylabel('Flow')
title('Select baseline time')
h=imrect;
pos = wait(h);
[~,bStart]=min(abs(time-pos(1)));
[~,bEnd]=min(abs(time-pos(1)-pos(3)));

title('Select response time')
h=imrect;
pos = wait(h);
[~,rStart]=min(abs(time-pos(1)));
[~,rEnd]=min(abs(time-pos(1)-pos(3)));

title('Selected baseline and response times')


figure
baseflowIndex=flowIndex(:,:,bStart:bEnd); 
mBaseflowIndex=squeeze(mean(baseflowIndex,3));

respflowIndex=flowIndex(:,:,rStart:rEnd)./mBaseflowIndex;
mRespflowIndex=squeeze(mean(respflowIndex,3));

% show flow map and select rectangular ROI
subplot(1,3,1)
imagesc(mBaseflowIndex)
caxis([prctile(mBaseflowIndex(:),5),prctile(mBaseflowIndex(:),95)]);
title('Baseline flow');

% show response map
subplot(1,3,2)
imagesc(mRespflowIndex)
caxis([prctile(mRespflowIndex(:),5),prctile(mRespflowIndex(:),95)]);
title('Response map, select ROI');
colorbar

for i=1:3
subplot(1,3,2)
h=imrect;
pos = wait(h);

pos=round(pos);
flowIndexROI=flowIndex(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
respROI=flowIndexROI./baseflowIndex(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3));
mRespRoi=squeeze(mean(mean(respROI,1),2));

%plot response normalized by baseline for the selected ROI
subplot(1,3,3)
hold on
plot(time,mRespRoi.*100)
hold off
xlabel('Time, s')
ylabel('Response, %')
title('ROI normalized flow')
end
subplot(1,3,2)
title('Response map');