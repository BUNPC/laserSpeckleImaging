% The script below identifies the type of the vessel (vein or artery)
% using the cardiac pulsatility information measured with LSCI.
%
% See "Cardiac pulsatility mapping and vessel type identification using
% laser speckle contrast imaging" by DD Postnov,S Erdener,K Kilic, and DA
% Boas, Biomedical Optics Express, 2018.
%
% This is only an example of the approach and can/should be improved to
% better utilize the pulsatility and BFI/contrast information.
%
% Additional m and MAT files are only used to load the raw data, calculate
% contrast and calculate powerSpectrums. They can be found at
% https://github.com/BUNPC/laserSpeckleImaging . They do not affect the
% approach demonstrated in the present code and can be replaced with
% corresponding functions/code.
%
% Other m-files required: readRLS, getSLSCI, getFFT, customColormaps
%
% MAT-files required: cmap
%
% Assumes section by section execution.
%
% Author: DD Postnov, PhD
% BOAS lab, Boston University
% BMI, Copenhagen University
% email address: dpostnov@sund.ku.dk / dpostnov@bu.edu
% Last revision: 1-November-2018
%%

%set directory with the current file and LSCI code folders
workdir='C:\Users\Dmitry\Documents\GitHub\laserSpeckleImaging';
cd(workdir);

addpath('Misc'); %for visualization
addpath('RawData'); %for reading and processing (LSCI) raw data
addpath('PowerAnalysis'); 

%load the custom colormaps
cmaps=customColormaps();

%set the file name
fName='E:\cardiac\2018614_crystal_aneth_1400x896x114.rls'; 

%Fourier processing parameters
pointsN=1024*4;
pointStart=100;
fftN=128;
harmonicsN=2;
averN=floor(pointsN./fftN);
mainFrqRange=[5,15];

%Vessel type identification parameters
%pP defines pixel belonging to parenchyma based on the
%power./median power relation. pP=1.2 seems to work well, arteries and
%veins normally have it around 5-6.
pP=1.2;
bN=100; %defines border of the image ignored in analysis


%Get data dimensions and sampling time to initialize arrays
[data,sampling,timeStamps]=readRLS(fName,0,1,[],'frame');
sampling=single(sampling)./1000;

meanBFI=zeros(size(data,1),size(data,2));
setH=zeros(size(data,1),size(data,2),averN,harmonicsN);
setP=zeros(size(data,1),size(data,2),averN,harmonicsN);
mainFrq=zeros(1,averN);


for i=1:1:averN
    waitbar(i./averN)
    
    %Contrast calculation from .rls file, replace with the appropriate
    %functions
    data=readRLS(fName,pointStart+(i-1)*fftN+1,fftN,[],'frame');
    data=getSLSCI(data,5,'gpu','none');
    maskNaN=isnan(data) | data<=0 | data>1; %find "bad" pixels
    data=1./(data.^2); % apply 1/K^2 model to get the BFI
    data(maskNaN==1)=0; % set "bad" pixels to zero flow
    
    %Get average BFI
    meanBFI=(meanBFI*(i-1)+squeeze(mean(data,3)))./i;
    
    %Get power spectrum, mean substracted
    [fftPow,~,f]=getFFT(data,sampling,fftN,'gpu');
    ts=squeeze(nanmean(nanmean(fftPow,1),2));
    
    %Identify the dominating frequency and set harmonics
    [~, idx]=max(ts);
    mainFrq(i)=f(idx);
    fListIdx=(1:1:harmonicsN).*idx;
    
    %If dominating frequency is within set cardiac range - get peak height
    %(H) and pedestal (P). Otherwise skip.
    if f(idx)>mainFrqRange(1) && f(idx)<mainFrqRange(2)
        for ii=1:1:length(fListIdx)
            powCardiac=fftPow(:,:,fListIdx(ii)-ceil((fListIdx(1)/3)*2):...
                fListIdx(ii)+floor((fListIdx(1)/3)*2));
            ts=squeeze(mean(mean(powCardiac,1),2));
            [~ , idx]=max(ts);
            setH(:,:,i,ii)=powCardiac(:,:,idx);
            setP(:,:,i,ii)=mean(powCardiac(:,:,[1:1:idx-2,idx+2:1:end]),3);
        end
    else
        disp('Cardiac is not dominant - strong motion is likely, skipping')
    end
end

%Example of approach to vessel identification using H and R information as
%shown in Fig 4, of "Cardiac pulsatility mapping and vessel type
%identification using laser speckle contrast imaging" article.

if sum(mainFrq)>0
    % Remove the data where cardiac activity was not dominant
    setP(:,:,toSkip==1,:)=[];
    setH(:,:,toSkip==1,:)=[];
    
    % Average the data, get relative and absolute peak hights
    P=squeeze(mean(setP,3));
    H=squeeze(mean(setH,3));
    R=H./P;
    H=H-P;
    
    typeMask=zeros(size(R));
    
    for i=1:1:harmonicsN
        % Process the first harmonic data
        tmpH=H(:,:,i);
        tmpR=R(:,:,i);
        tmpP=P(:,:,i);
        
        % Get median values ~ parenchyma
        medH=imresize(tmpH,0.1);
        medH=medfilt2(medH,[21 21]);
        medH=imresize(medH,size(tmpH));
        medR=imresize(tmpR,0.1);
        medR=medfilt2(medR,[21 21]);
        medR=imresize(medR,size(tmpH));
        medP=imresize(tmpP,0.1);
        medP=medfilt2(medP,[21 21]);
        medP=imresize(medP,size(tmpP));
        
        mapH=tmpH./medH;
        mapR=tmpR./medR;
        mapP=tmpP./medP;
        
        mapP=imgaussfilt(mapP,1);
        mapH=imgaussfilt(mapH,1);
        mapR=imgaussfilt(mapR,1);
        
        
        
        bmask=zeros(size(tmpP));
        bmask(bN+1:1:end-bN,bN+1:1:end-bN)=1;
        
        %define the deviation of R and P for parenchyma
        rP=nanstd(mapR(mapP<=pP & bmask==1))*2;
        hP=nanstd(mapH(mapP<=pP & bmask==1))*2;
        
        %Find likely values for of R and H for arteries and veins. Logic
        %is following:
        %For veins R and H are likely to be lower than of parenchyma or
        %arteries, while the power pedestal is much higher than for
        %parenchyma (if BFI is used for power spectrum calculations).
        %For arteries it is similar but both R and H are likely to be
        %higher than those of parenchyma or veins.
        %WORKS ONLY IF BOTH TYPES ARE PRESENT
        rA=prctile(mapR(bmask==1 & mapP>pP),90);
        rV=prctile(mapR(bmask==1 & mapP>pP),10);
        hA=prctile(mapH(bmask==1 & mapP>pP),90);
        hV=prctile(mapH(bmask==1 & mapP>pP),10);
        
        %If only one type of the vessel is present - following logic is
        %recommended: mapR>0 & mapH>0 - likely artery; mapR>0 - maybe
        %artery; mapR<0 - maybe vein; mapR<0 & mapH<0 - likely vein
        
        %Separation value for R
        rS=1.0*(rA+rV)./2;
        
        map=(mapR-rS)./(rA-rS);
        map(map>1)=1; %maybe artery
        map(map<-1)=-1; %maybe vein
        map(mapR>=1+rP)=2; %likely artery
        map(mapH<=1-hP)=-2; %likely vein
        map(mapP<=pP)=-4; %likely parenchyma
        typeMask(:,:,i)=map;
    end
    
    %Combine results from harmonics
    typeMask=squeeze(mean(typeMask,3));
    typeMask(typeMask>1 & typeMask<2)=2;
    typeMask(typeMask<-1 & typeMask>-2)=-2;
    
    %Display
    figure
    imagesc(meanBFI);
    caxis([prctile(meanBFI(:),1),prctile(meanBFI(:),99)])
    colormap(cmaps.imgColor)
    title('mean BFI')
    
    figure
    imagesc(typeMask);
    caxis([-4,4])
    colormap(cmaps.twoType)
    title('vessel type')
end

