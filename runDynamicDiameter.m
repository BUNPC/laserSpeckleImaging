%Example of a run script for the diameter and flow estimation from the 
%blood flow imaging data (i.e. LSCI).
%
% For details see
% Postnov DD, Tuchin VV, Sosnovtseva O. Estimation of vessel diameter and 
% blood flow dynamics from laser speckle images. Biomedical optics express.
% 2016 Jul 1;7(7):2759-68.
%
% Assumes section by section execution.
%
% Author: DD Postnov, PhD
% BOAS lab, Boston University
% BMI, Copenhagen University
% email address: dpostnov@sund.ku.dk / dpostnov@bu.edu
% Last revision: 1-November-2018

%% Add path
%set directory with the current file and LSCI code folders
workdir='C:\Users\Dmitry\Documents\GitHub\laserSpeckleImaging';
cd(workdir);

addpath('Misc'); %for visualization
addpath('Raw'); %for reading and processing (LSCI) raw data
addpath('Segmentation'); 

%load the custom colormaps
cmaps=customColormaps();

%set name and folder of the raw data file
fileName='20180903__109WT_male_Leftwhiskerstimulation_100to300';
folderName='E:\Dropbox\Upload raw files\';
fileNameFull=[folderName,fileName,'.rls'];

%short name to be used in Results structure
expName='prop';

disp('Path loaded');

%% Read and process the .rls raw data.
% Replace with appropriate functions if the raw data is not in .rls format
data=readRLS(fileNameFull,0,0,[200,800;200,800],'rect');
% Do the temporal contrast analysis/replace with getSLSCI if appropriate
LSCI=getTLSCI(data,25,'gpu','kernel');

disp('Data processed');

%% Prepare the segmentation data 

dataSegm=LSCI;
% Clean up extremely dark/bright regions
imgMean=mean(dataSegm,3);
dataSegm(dataSegm<prctile(imgMean(:),0.1))=prctile(imgMean(:),0.1);
dataSegm(dataSegm>prctile(imgMean(:),99.9))=prctile(imgMean(:),99.9);

% Make foreground bright and convert to the flow imaging (for the contrast data data)
dataSegm=1./(dataSegm.*dataSegm);

% Set data dimensions
sizeT=size(dataSegm,3);
sizeX=size(dataSegm,1);
sizeY=size(dataSegm,2);

disp('Segmentation data prepared');


%% Do the segmentation
% Select vessel ROIs
img=squeeze(mean(dataSegm,3));
figure;
displayImg(img,[1,99],'Select ROIs',cmaps.imgColor);
ROIs=getTrapezoidSelection(img);

%Segmentation parameters
expD=5; %Exepcted minimal diameter. Affects various spatial clean ups.
smoothS1=1000; % N frames for temporal filtering for secondary guess. Does not affect temporal resolution. No periodic vessel shift should happen during this time.
smoothS2=1000; % N frames for the cleanup filter of the secondary guess. Same as smoothS1.
smoothS3=100; % N fames for temporal filtering raw data. DOES affect temporal resolution. Any dynamics faster than this number won't be resolvable.
closeLoops1=1; % N of close operation loops (dilation+erosion) before the largest object selection. Always small - default is 1.
closeLoops2=3; % N of close operation loops (dilation+erosion) after the largest object selection. Cleans up filling artifacts. Large values won't harm straight vessels but may affect curved vessels. Default is 3.
medfiltS1=0; %N pixels for directional median filtering of the raw data. Affects the spatial resolution (decreases shaprness). Should be applied when extremely dark/bright artifacts are present. Default is 0.
medfiltS2=expD; %N pixels in median filtering of the centerline location. Not applied when shape='straight'. Larger values reduce artifacts but will be harmfull for strongly curved vessels. Default is expD.
medfiltS3=5; %N pixels for median filtering of the profile prior to averaging. Reduces noise artifacts. Larger values might be harmfull for strongly curved vessels. Default value is 3 or 5.
medfiltS4=100; %N frames for filtering of the estimated diameter. DOES affect the temporal resolution. Any dynamics faster than this number won't be resolvable.
shape='curved'; % 'straight' when vessels are surely straight, otherwise anything else, i.e. 'curved'.
method='mvar'; %Diameter estimation method, only 'mvar' is availible atm.


%Allocating space for average mask over the whole region
Result.(sprintf(char(expName))).segmMaskMean=zeros(sizeX,sizeY);
Result.(sprintf(char(expName))).segmMaskROIs=zeros(sizeX,sizeY);
Result.(sprintf(char(expName))).lineMaskMean=zeros(sizeX,sizeY);
%Segmentation in selected ROI's
for i=1:1:length(ROIs)
    tic
    disp(['Processing ROI ',num2str(i),' - started']);
    x=ROIs(i).x;
    y=ROIs(i).y;
    mask=ROIs(i).mask;
    dataROI=zeros(x(2)-x(1)+1,y(2)-y(1)+1,sizeT,'single');
    for t=1:1:sizeT
        dataROI(:,:,t)=dataSegm(x(1):x(2),y(1):y(2),t).*mask+squeeze(mean(mean(dataSegm(x(1):x(2),y(1):y(2),t),1),2)).*(single(ones(x(2)-x(1)+1,y(2)-y(1)+1))-mask);
    end
    region=twoPhaseSegmentation(dataROI,expD,smoothS1,smoothS2,smoothS3,closeLoops1,closeLoops2,medfiltS1,medfiltS2,medfiltS3,medfiltS4,shape,method);
    
    %Results structure
    Result.(sprintf(char(expName))).(sprintf('r%d',i))= region;
    Result.(sprintf(char(expName))).segmMaskMean(x(1):x(2),y(1):y(2))=Result.(sprintf(char(expName))).segmMaskMean(x(1):x(2),y(1):y(2))+squeeze(mean(region.maskVessel3,3));
    Result.(sprintf(char(expName))).lineMaskMean(x(1):x(2),y(1):y(2))=Result.(sprintf(char(expName))).lineMaskMean(x(1):x(2),y(1):y(2))+squeeze(mean(region.maskCenterLine,3));
    Result.(sprintf(char(expName))).segmMaskROIs(x(1):x(2),y(1):y(2))=Result.(sprintf(char(expName))).segmMaskROIs(x(1):x(2),y(1):y(2))+(squeeze(mean(region.maskVessel3,3))>0.5).*i;
    
    %Get flow and speed calculations
    dataROI=LSCI(x(1):x(2),y(1):y(2),:);
    dataROI=1./(dataROI.*dataROI);
    if sizeT>smoothS3 && smoothS3>1
    for x=1:1:size(dataROI,1)
        for y=1:1:size(dataROI,2)
            dataROI(x,y,:)=smooth(dataROI(x,y,:),smoothS3);
        end
    end
    end 
    Result.(sprintf(char(expName))).(sprintf('r%d',i)).speedLine=squeeze(nansum(nansum(region.maskCenterLine.*dataROI,1),2)./sum(sum(region.maskCenterLine,1),2));
    Result.(sprintf(char(expName))).(sprintf('r%d',i)).speedVessel=squeeze(nansum(nansum(region.maskVessel3.*dataROI,1),2)./sum(sum(region.maskVessel3,1),2));
    Result.(sprintf(char(expName))).(sprintf('r%d',i)).flowLine=squeeze(nansum(nansum(region.maskCenterLine.*dataROI,1),2)./sum(sum(region.maskCenterLine,1),2)).*(region.d.*region.d)';
    Result.(sprintf(char(expName))).(sprintf('r%d',i)).flowVessel=squeeze(nansum(nansum(region.maskVessel3.*dataROI,1),2)./sum(sum(region.maskVessel3,1),2)).*(region.d.*region.d)';    
    
    %Save the dataROI if required
    %Result.(sprintf(char(expName))).(sprintf('r%d',i)).dataROI=dataROI;
    clearvars region;
    disp(['Processing ROI ',num2str(i),' - finished']);
    toc
end

%% Example of the visualisation:
figure
subplot (1,2,1)
imgMean=1./LSCI;
imgMean=squeeze(mean(imgMean.*imgMean,3));
displayImg(imgMean,[1,99],'average BFI',cmaps.imgColor);
hold on
visboundaries(mean(Result.(sprintf(char(expName))).segmMaskMean,3)>0.5)
hold off
title('BFI map and segmentation')
subplot (1,2,2)
imagesc(Result.(sprintf(char(expName))).segmMaskROIs);
title('Segmentation ROIs')

windows=ceil(length(ROIs)./3);
panels=3;
for i=1:1:windows
    figure
    for j=1:1:min(panels,length(ROIs)-(i-1)*panels)
        subplot(3,panels,j)
        plot(1:1:sizeT,Result.(sprintf(char(expName))).(sprintf('r%d',(i-1)*panels+j)).d)
        title(['ROI ',num2str((i-1)*panels+j),' d, pixels'])
        subplot(3,panels,j+panels)
        plot(1:1:sizeT,Result.(sprintf(char(expName))).(sprintf('r%d',(i-1)*panels+j)).speedVessel)
        title(['ROI ',num2str((i-1)*panels+j),' speed, AU'])
        subplot(3,panels,j+2*panels)
        plot(1:1:sizeT,Result.(sprintf(char(expName))).(sprintf('r%d',(i-1)*panels+j)).flowVessel)
        title(['ROI ',num2str((i-1)*panels+j),' flow, AU'])
    end
end