%Example of the run script for the laser speckle contrast analysis.
%Includes loading the raw data, contrast calculation, some visualization
%
% Other m-files required: getLSCIFolderIMG, getSLSCI, displayImg,
% customColormaps
%
% MAT-files required: cmaps
%
% Assumes section by section execution.
%
% Author: DD Postnov, PhD
% BOAS lab, Boston University
% BMI, Copenhagen University
% email address: dpostnov@sund.ku.dk / dpostnov@bu.edu
% Last revision: 1-November-2018

%% Read and process at the same time - the most standard procedure

%set directory with the current file and LSCI code folders
workdir='C:\Users\Dmitry\Documents\GitHub\laserSpeckleImaging';
cd(workdir);

addpath('Misc'); %for visualization
addpath('RawData'); %for reading and processing (LSCI) raw data

folderName='D:\Evren\trial 1 left'; % folder with raw files
burstStart=1; % first burst to read
burstN=[]; % number of bursts to read, [] for all bursts
kernelSize=7; % kernel size for LSCI processing. 7 or 5 for sLSCI, 25 for tLSCI
lscType='sLSCI'; % type of contrast analysis
procTypeLSCI='gpu'; %processing type for contrast analysis
dsType='none'; %data downsampling (none or kernel)
procType='cluster'; % single thread or parallel processing of multiple raw files

%Load and process the raw data recorded with A.Dunn software (not .rls)
[LSCI,time,imgsPerBurst]=getLSCIFolderIMG(procType,folderName,burstStart,burstN,lscType,kernelSize,procTypeLSCI,dsType);
% substract first time value (in case if no loop over burst files is used)
time=time-time(1);
% save results
save([folderName,'\LSCI','.mat'],'LSCI','time','lscType','dsType','imgsPerBurst','-v7.3');

%% Convert contrast to flow
close all
% assuming that flow is 1/K^2
flowIndex=1./(LSCI.*LSCI);
figure
subplot(2,2,2)
plot(time,squeeze(mean(mean(LSCI,1),2)));
xlabel('Time,s')
ylabel('Contrast')
title('Select baseline time')
h=imrect;
pos = wait(h);
[~,bStart]=min(abs(time-pos(1)));
[~,bEnd]=min(abs(time-pos(1)-pos(3)));
if bStart<1
    bStart=1;
end
if bEnd>size(flowIndex,3)
    bEnd=size(flowIndex,3);
end

title('Select response time')
h=imrect;
pos = wait(h);
[~,rStart]=min(abs(time-pos(1)));
[~,rEnd]=min(abs(time-pos(1)-pos(3)));
if rStart<1
    rStart=1;
end
if rEnd>size(flowIndex,3)
    rEnd=size(flowIndex,3);
end
title('Selected baseline and response times')

baseflowIndex=flowIndex(:,:,bStart:bEnd); 
mBaseflowIndex=squeeze(mean(baseflowIndex,3));

respflowIndex=flowIndex(:,:,rStart:rEnd)./mBaseflowIndex;
mRespflowIndex=squeeze(mean(respflowIndex,3));

% show flow map and select rectangular ROI
subplot(2,2,1)
displayImg(mBaseflowInde,[1,99],'baseline BFI',cmaps.imgColor);

% show response map
subplot(2,2,3)
displayImg(mRespflowIndex,[1,99],'Response map, select ROI',cmaps.imgColor);

%select and plot ROIs
for i=1:3
subplot(2,2,3)
h=imrect;
pos = wait(h);
delete(h);
hold on
plot([pos(1),pos(1),pos(1)+pos(3),pos(1)+pos(3),pos(1)],[pos(2),pos(2)+pos(4),pos(2)+pos(4),pos(2),pos(2)],'LineWidth',3);
hold off

pos=round(pos);
flowIndexROI=flowIndex(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
respROI=flowIndexROI./baseflowIndex(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3));
mRespRoi=squeeze(mean(mean(respROI,1),2));

subplot(2,2,4)
hold on
plot(time,mRespRoi.*100)
hold off
xlabel('Time, s')
ylabel('Response, %')
title('ROI normalized flow')
end
subplot(2,2,3)
title('Response map');