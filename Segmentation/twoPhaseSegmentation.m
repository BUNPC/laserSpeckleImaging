%twoPhaseSegmentation - identifies straight or slightly curved tube-like 
%object (vessel) in the 3D (X,Y,time) data with bright foreground.
%
% Syntax:  region=twoPhaseSegmentation(data,expD,smoothS1,smoothS2,
%                   smoothS3,closeLoops1,closeLoops2,medfiltS1,medfiltS2,
%                   medfiltS3,medfiltS4,shape,method)
%
% Inputs:
%    data        - 2D (XY) or 3D (XYT) data
%    expD        - Exepcted minimal diameter. Affects spatial cleanups.
%    smoothS1    - N frames for temporal filtering for secondary guess. 
%                  Does not affect temporal resolution. No periodic vessel
%                  shift should happen during this time.
%    smoothS2    - N frames for the cleanup filter of the secondary guess.
%                  Same as smoothS1.
%    smoothS3    - N fames for temporal filtering raw data. 
%                  DOES affect temporal resolution. Any dynamics faster 
%                  than this number won't be resolvable.
%    closeLoops1 - N of close operation loops (dilation+erosion) before the
%                  largest object selection. Recommended value is 1.
%    closeLoops2 - N of close operation loops (dilation+erosion) after the
%                  largest object selection. Cleans up filling artifacts.
%                  Large values won't harm straight vessels but may affect
%                  curved vessels. Default is 3.
%    medfiltS1   - N pixels for directional median filtering of the raw
%                  data. Reduces spatial resolution. Should be applied when
%                  extremely dark/bright artifacts are present. Recommended
%                  value is 1.
%    medfiltS2   - N pixels in median filtering of the centerline location.
%                  Not applied when shape='straight'. Larger values reduce
%                  artifacts but will be harmfull for strongly curved
%                  vessels. Recommended value is the same as expD.
%    medfiltS3   - N pixels for median filtering of the profile prior to
%                  averaging. Reduces noise artifacts. Larger values might
%                  be harmfull for strongly curved vessels. Recommended
%                  value is 3 or 5.
%    medfiltS4   - N frames for filtering of the estimated diameter. DOES 
%                  affect the temporal resolution. Any dynamics faster than
%                  this number won't be resolvable.
%    shape       - 'straight' when vessels are surely straight, otherwise
%                  anything else, i.e. 'curved'.
%    method      - Object boundaries detection method, 
%                  only 'mvar' is availible atm.
%
% Outputs:
%    region      - structure that contains segmentation data and diameter
%                  measurements
%
% Example:
%    region=twoPhaseSegmentation(data,5,50,50,10,1,3,0,5,3,3,'straight'
%               ,'mvar');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: DD Postnov, PhD
% BOAS lab, Boston University
% BMI, Copenhagen University
% email address: dpostnov@sund.ku.dk / dpostnov@bu.edu
% Last revision: 1-November-2018

%------------- BEGIN CODE --------------

function region=twoPhaseSegmentation(data,expD,smoothS1,smoothS2,smoothS3,closeLoops1,closeLoops2,medfiltS1,medfiltS2,medfiltS3,medfiltS4,shape,method)

%% Prepare data in the region and meta data
data=single(data);
sizeY=size(data,1);
sizeX=size(data,2);
sizeT=size(data,3);

roughnessLimit=0.5;
filter1SD=expD/6;
filter2SD=expD/12;

region.input.expD=expD;
region.input.smoothS1=smoothS1;
region.input.smoothS2=smoothS2;
region.input.smoothS3=smoothS3;
region.input.closeLoops1=closeLoops1;
region.input.closeLoops2=closeLoops2;
region.input.medfiltS1=medfiltS1;
region.input.medfiltS2=medfiltS2;
region.input.medfiltS3=medfiltS3;
region.input.medfiltS4=medfiltS4;
region.input.shape=shape;
region.input.method=method;

%%
dataS=zeros(sizeY,sizeX,sizeT,'single');
padval=squeeze(nanmean(nanmean(data,1),2))';
for t=1:1:sizeT
    dataS(:,:,t)=imgaussfilt((squeeze(data(:,:,t))),filter1SD);
end
meanImage=squeeze(mean(dataS(:,:,:),3));
if sizeT>smoothS1 && smoothS1>1
    for x=1:1:sizeX
        for y=1:1:sizeY
            dataS(y,x,:)=smooth(dataS(y,x,:),smoothS1);
        end
    end
end

%% Allocate masks
maskVesselV=zeros(sizeY,sizeX,sizeT,'logical');
maskVesselH=zeros(sizeY,sizeX,sizeT,'logical');
maskLine=zeros(sizeY,sizeX,sizeT,'logical');

%% Allocate some of results
DleftH=zeros(sizeT,sizeY);
DrightH=zeros(sizeT,sizeY);
DleftV=zeros(sizeT,sizeX);
DrightV=zeros(sizeT,sizeX);
alpha=zeros(1,sizeT);

%% Run check in horizontal direction (vessel is vertical)
for y=1:1:sizeY
    %Initial guess based on the mean data
    yy=meanImage(y,:);
    [~,imax]=max(yy);
    
    icur=imax;
    dif=1;
    while dif>0 && icur>1
        dif=yy(icur)-yy(icur-1);
        icur=icur-1;
    end
    yyleft=yy(icur);
    
    icur=imax;
    dif=1;
    while dif>0 && icur<length(yy)
        dif=yy(icur)-yy(icur+1);
        icur=icur+1;
    end
    yyright=yy(icur);
    
    level=max(yyleft,yyright);
    icur=imax;
    while yy(icur)>level
        icur=icur-1;
    end
    ileftGuess=icur;
    
    icur=imax;
    while yy(icur)>level
        icur=icur+1;
    end
    irightGuess=icur;
    
    % start each frame processing
    for t=1:1:sizeT
        yy=squeeze(dataS(y,:,t));
        yyMidLevel=mean(yy);
        yyGuess=yy(ileftGuess:irightGuess);
        [~,imax]=max(yyGuess);
        imax=imax+ileftGuess-1;
        
        icur=imax;
        dif=1;
        while dif>0 && icur>1
            if yy(icur-1)<yyMidLevel
                dif=yy(icur)-yy(icur-1);
            end
            icur=icur-1;
        end
        ileftMax=icur;
        
        icur=imax;
        dif=1;
        while dif>0 && icur<length(yy)
            if yy(icur+1)<yyMidLevel
                dif=yy(icur)-yy(icur+1);
            end
            icur=icur+1;
        end
        irightMin=icur;
        
        level=max(yy(ileftMax),yy(irightMin));
        icur=imax;
        while yy(icur)>level && icur>1
            icur=icur-1;
        end
        ileft=icur;
        
        icur=imax;
        while yy(icur)>level && icur<length(yy)
            icur=icur+1;
        end
        iright=icur;
        
        maskVesselH(y,ileft:iright,t)=1;
        DleftH(t,y)=ileft;
        DrightH(t,y)=iright;
    end
end

region.maskVessel1H=maskVesselH;

%% Run check in vertical direction (vessel is horizontal)

for x=1:1:sizeX
    
    %Initial guess based on the mean data
    xx=meanImage(:,x);
    [~,imax]=max(xx);
    
    icur=imax;
    dif=1;
    while dif>0 && icur>1
        dif=xx(icur)-xx(icur-1);
        icur=icur-1;
    end
    xxleft=xx(icur);
    
    icur=imax;
    dif=1;
    while dif>0 && icur<length(xx)
        dif=xx(icur)-xx(icur+1);
        icur=icur+1;
    end
    xxright=xx(icur);
    
    level=max(xxleft,xxright);
    icur=imax;
    while xx(icur)>level
        icur=icur-1;
    end
    ileftGuess=icur;
    
    icur=imax;
    while xx(icur)>level
        icur=icur+1;
    end
    irightGuess=icur;
    
    % start each frame processing
    for t=1:1:sizeT
        xx=squeeze(dataS(:,x,t));
        xxMidLevel=mean(xx);
        xxGuess=xx(ileftGuess:irightGuess);
        [~,imax]=max(xxGuess);
        imax=imax+ileftGuess-1;
        
        icur=imax;
        dif=1;
        while dif>0 && icur>1
            if xx(icur-1)<xxMidLevel
                dif=xx(icur)-xx(icur-1);
            end
            icur=icur-1;
        end
        ileftMax=icur;
        
        icur=imax;
        dif=1;
        while dif>0 && icur<length(xx)
            if xx(icur+1)<xxMidLevel
                dif=xx(icur)-xx(icur+1);
            end
            icur=icur+1;
        end
        irightMin=icur;
        
        level=max(xx(ileftMax),xx(irightMin));
        icur=imax;
        while xx(icur)>level && icur>1
            icur=icur-1;
        end
        ileft=icur;
        
        icur=imax;
        while xx(icur)>level && icur<length(xx)
            icur=icur+1;
        end
        iright=icur;
        
        maskVesselV(ileft:iright,x,t)=1;
        DleftV(t,x)=ileft;
        DrightV(t,x)=iright;
        
    end
    
end
region.maskVessel1V=maskVesselV;

%% Feel in the "time holes" in both directions.
maskTMP1=zeros(size(maskVesselH),'single');
maskTMP2=zeros(size(maskVesselV),'single');

for t=1:1:sizeT
    tLeft=t-floor(smoothS2/2);
    tRight=t+floor(smoothS2/2);
    if tLeft<1
        tLeft=1;
    end
    if tRight>sizeT
        tRight=sizeT;
    end
    maskTMP1(:,:,t)=squeeze(mean(single(maskVesselH(:,:,tLeft:tRight)),3));
    maskTMP2(:,:,t)=squeeze(mean(single(maskVesselV(:,:,tLeft:tRight)),3));
    
end
maskVesselH=maskTMP1>0.5;
maskVesselV=maskTMP2>0.5;

%% Clean up small roughness in both directions

maskTMP=single(maskVesselH);
for t=1:1:sizeT
    maskTMP(:,:,t)=imgaussfilt((squeeze(maskTMP(:,:,t))),filter2SD);
end
maskVesselH=maskTMP>roughnessLimit;

maskTMP=single(maskVesselV);
for t=1:1:sizeT
    maskTMP(:,:,t)=imgaussfilt((squeeze(maskTMP(:,:,t))),filter2SD);
end
maskVesselV=maskTMP>roughnessLimit;

%%
SE=[0 1 0; 1 1 1; 0 1 0];
if closeLoops1>0
    for t=1:1:sizeT
        
        %Horizontal
        maskTMP=padarray(squeeze(maskVesselH(:,:,t)),[closeLoops1+1,closeLoops1+1],0);
        for i=1:1:closeLoops1
            maskTMP=imdilate(maskTMP,SE) ;
        end
        for i=1:1:closeLoops1
            maskTMP=imerode(maskTMP,SE) ;
        end
        maskVesselH(:,:,t)=maskTMP(closeLoops1+2:end-closeLoops1-1,closeLoops1+2:end-closeLoops1-1);
        
        %Vertical
        maskTMP=padarray(squeeze(maskVesselV(:,:,t)),[closeLoops1+1,closeLoops1+1],0);
        for i=1:1:closeLoops1
            maskTMP=imdilate(maskTMP,SE) ;
        end
        for i=1:1:closeLoops1
            maskTMP=imerode(maskTMP,SE) ;
        end
        maskVesselV(:,:,t)=maskTMP(closeLoops1+2:end-closeLoops1-1,closeLoops1+2:end-closeLoops1-1);
    end
end

%% Deterimine the largest object and remove all others in both directions, choose direction
for t=1:1:sizeT
    %Horizontal
    CC = bwconncomp(squeeze(maskVesselH(:,:,t)), 4);
    S = regionprops(CC, 'Area');
    L = labelmatrix(CC);
    [~,idx]=max([S.Area]);
    maskVesselH(:,:,t) = ismember(L, idx);
    
    %Vertical
    CC = bwconncomp(squeeze(maskVesselV(:,:,t)), 4);
    S = regionprops(CC, 'Area');
    L = labelmatrix(CC);
    [~,idx]=max([S.Area]);
    maskVesselV(:,:,t) = ismember(L, idx);
end

%% Do erosion/dilation loops over all frames
SE=[0 1 0; 1 1 1; 0 1 0];
if closeLoops2>0
    for t=1:1:sizeT
        
        %Horizontal
        maskTMP=padarray(squeeze(maskVesselH(:,:,t)),[closeLoops2+1,closeLoops2+1],0);
        for i=1:1:closeLoops2
            maskTMP=imdilate(maskTMP,SE) ;
        end
        for i=1:1:closeLoops2
            maskTMP=imerode(maskTMP,SE) ;
        end
        maskVesselH(:,:,t)=maskTMP(closeLoops2+2:end-closeLoops2-1,closeLoops2+2:end-closeLoops2-1);
        
        %Vertical
        maskTMP=padarray(squeeze(maskVesselV(:,:,t)),[closeLoops2+1,closeLoops2+1],0);
        for i=1:1:closeLoops2
            maskTMP=imdilate(maskTMP,SE) ;
        end
        for i=1:1:closeLoops2
            maskTMP=imerode(maskTMP,SE) ;
        end
        maskVesselV(:,:,t)=maskTMP(closeLoops2+2:end-closeLoops2-1,closeLoops2+2:end-closeLoops2-1);
    end
end

%%
% Choose directions with larger (on average) object in it
if (sum(maskVesselH(:))>sum(maskVesselV(:)))
    chosenDirection='H';
    region.maskVessel2=maskVesselH;
else
    chosenDirection='V';
    region.maskVessel2=maskVesselV;
end
region.sDir=chosenDirection;

%% Smooth the raw data
if sizeT>medfiltS1 && medfiltS1>1
    if (strcmp(chosenDirection,'H'))
        for t=1:1:sizeT
            for y=1:1:sizeY
                data(y,:,t)=medfilt1(squeeze(data(y,:,t)),medfiltS1,'omitnan','truncate');
            end
        end
    else
        for t=1:1:sizeT
            for x=1:1:sizeX
                data(:,x,t)=medfilt1(squeeze(data(:,x,t)),medfiltS1,'omitnan','truncate');
            end
        end
    end
end

if sizeT>smoothS3 && smoothS3>1
    for x=1:1:sizeX
        for y=1:1:sizeY
            data(y,x,:)=smooth(data(y,x,:),smoothS3);
        end
    end
end

%% Initial centerline detection
setTS1=cell(1,sizeT);
setTS2=cell(1,sizeT);
nsse=zeros(1,sizeT);
if (strcmp(chosenDirection,'H'))
    for t=1:1:sizeT
        k=1;
        indTS1=zeros(1,sum(squeeze(sum(squeeze(maskVesselH(:,:,t)),2))>0));
        indTS2=zeros(1,sum(squeeze(sum(squeeze(maskVesselH(:,:,t)),2))>0));
        width=zeros(1,sum(squeeze(sum(squeeze(maskVesselH(:,:,t)),2))>0));
        for y=1:1:sizeY
            ts=single(squeeze(maskVesselH(y,:,t)));
            sumTS=sum(ts);
            if sumTS>0
                indTS1(k)=y;
                indTS2(k)=find(ts==1,1)+floor(sum(ts)/2);
                width(k)=sum(ts);
                k=k+1;
            end
        end
        if strcmp(shape,'straight')
            p=polyfit(indTS1,indTS2,1);
            nsse(t)=mean(((p(1).*indTS1+p(2) - indTS2).^2));
            indTS1=1:1:sizeY;
            indTS2=p(1).*indTS1+p(2);
        elseif medfiltS2>1
            indTS2 = medfilt1(indTS2,medfiltS2,'omitnan','truncate');
        end
        indTS2=round(indTS2);
        indTS2(indTS2<1)=1;
        indTS2(indTS2>sizeX)=sizeX;
        idx1d=(indTS2-1).*sizeY+indTS1;
        tmpMask=zeros(sizeY,sizeX);
        tmpMask(idx1d)=1;
        S = regionprops(tmpMask, 'Orientation');
        
        maskLine(:,:,t)=tmpMask;
        alpha(t)=S.Orientation;
        setTS1{t}=indTS1;
        setTS2{t}=indTS2;
        
    end
else
    for t=1:1:sizeT
        k=1;
        indTS1=zeros(1,sum(squeeze(sum(squeeze(maskVesselV(:,:,t)),1))>0));
        indTS2=zeros(1,sum(squeeze(sum(squeeze(maskVesselV(:,:,t)),1))>0));
        width=zeros(1,sum(squeeze(sum(squeeze(maskVesselV(:,:,t)),1))>0));
        for x=1:1:sizeX
            ts=single(squeeze(maskVesselV(:,x,t)));
            sumTS=sum(ts);
            if sumTS>0
                indTS1(k)=x;
                indTS2(k)=find(ts==1,1)+floor(sum(ts)/2);
                width(k)=sum(ts);
                k=k+1;
            end
        end
        if strcmp(shape,'straight')
            p=polyfit(indTS1,indTS2,1);
            nsse(t)=mean(((p(1).*indTS1+p(2) - indTS2).^2));
            indTS1=1:1:sizeX;
            indTS2=p(1)*indTS1+p(2);
        elseif medfiltS2>1
            indTS2 = medfilt1(indTS2,medfiltS2,'omitnan','truncate');
        end
        indTS2=round(indTS2);
        indTS2(indTS2<1)=1;
        indTS2(indTS2>sizeY)=sizeY;
        idx1d=(indTS1-1).*sizeY+indTS2;
        tmpMask=zeros(sizeY,sizeX);
        tmpMask(idx1d)=1;
        S = regionprops(tmpMask, 'Orientation');
        
        maskLine(:,:,t)=tmpMask;
        alpha(t)=S.Orientation;
        setTS1{t}=indTS1;
        setTS2{t}=indTS2;
    end
end

region.nsse=nsse;
%% Calculate the profile

if (strcmp(chosenDirection,'H'))
    profile=ones(sizeX*2+1,sizeT).*padval;
    for t=1:1:sizeT
        indTS1=setTS1{t};
        indTS2=setTS2{t};
        vals=zeros(1,length(indTS1));
        for y=1:1:length(indTS1)
            vals(y)=data(indTS1(y),indTS2(y),t);
        end
        if medfiltS3>1
            vals=medfilt1(vals,medfiltS3,'omitnan','truncate');
        end
        profile(sizeX+1,t)=nanmean(vals);
        maxShiftPlus=min(sizeX-min(indTS2),sizeX);
        maxShiftMinus=min(max(indTS2),sizeX);
        for kk=1:1:maxShiftPlus
            indTS2tmp=indTS2+kk;
            indTS1tmp=indTS1;
            indTS1tmp(indTS2tmp>sizeX)=[];
            indTS2tmp(indTS2tmp>sizeX)=[];
            vals=zeros(1,length(indTS1tmp));
            for y=1:1:length(indTS1tmp)
                vals(y)=data(indTS1tmp(y),indTS2tmp(y),t);
            end
            if medfiltS3>1
                vals=medfilt1(vals,medfiltS3,'omitnan','truncate');
            end
            profile(sizeX+1+kk,t)=nanmean(vals);
        end
        
        for kk=1:1:maxShiftMinus
            indTS2tmp=indTS2-kk;
            indTS1tmp=indTS1;
            indTS1tmp(indTS2tmp<1)=[];
            indTS2tmp(indTS2tmp<1)=[];
            vals=zeros(1,length(indTS1tmp));
            for y=1:1:length(indTS1tmp)
                vals(y)=data(indTS1tmp(y),indTS2tmp(y),t);
            end
            if medfiltS3>1
                vals=medfilt1(vals,medfiltS3,'omitnan','truncate');
            end
            profile(sizeX+1-kk,t)=nanmean(vals);
        end
        ts=profile(:,t);
        ts(isnan(ts))=padval(t);
        profile(:,t)=ts;
    end
    
else
    profile=ones(sizeY*2+1,sizeT).*padval;
    for t=1:1:sizeT
        indTS1=setTS1{t};
        indTS2=setTS2{t};
        vals=zeros(1,length(indTS1));
        for y=1:1:length(indTS1)
            vals(y)=data(indTS2(y),indTS1(y),t);
        end
        if medfiltS3>1
            vals=medfilt1(vals,medfiltS3,'omitnan','truncate');
        end
        profile(sizeY+1,t)=nanmean(vals);
        maxShiftPlus=min(sizeY-min(indTS2),sizeY);
        maxShiftMinus=min(max(indTS2),sizeY);
        
        for kk=1:1:maxShiftPlus
            indTS2tmp=indTS2+kk;
            indTS1tmp=indTS1;
            indTS1tmp(indTS2tmp>sizeY)=[];
            indTS2tmp(indTS2tmp>sizeY)=[];
            vals=zeros(1,length(indTS1tmp));
            for y=1:1:length(indTS1tmp)
                vals(y)=data(indTS2tmp(y),indTS1tmp(y),t);
            end
            if medfiltS3>1
                vals=medfilt1(vals,medfiltS3,'omitnan','truncate');
            end
            profile(sizeY+1+kk,t)=nanmean(vals);
        end
        for kk=1:1:maxShiftMinus
            indTS2tmp=indTS2-kk;
            indTS1tmp=indTS1;
            indTS1tmp(indTS2tmp<1)=[];
            indTS2tmp(indTS2tmp<1)=[];vals=zeros(1,length(indTS1tmp));
            for y=1:1:length(indTS1tmp)
                vals(y)=data(indTS2tmp(y),indTS1tmp(y),t);
            end
            if medfiltS3>1
                vals=medfilt1(vals,medfiltS3,'omitnan','truncate');
            end
            profile(sizeY+1-kk,t)=nanmean(vals);
        end
        ts=profile(:,t);
        ts(isnan(ts))=padval(t);
        profile(:,t)=ts;
    end
end

%% Clean up profile

%% Get the vessel bordere based on FWHM

%% Get the vessel bordere based on second gradient zero

%% Get the vessel borderes based on the minimzed class variation
if strcmp(method,'mvar')
    interpStep=0.1;
    
    idxIni=zeros(1,sizeT);
    idxL=zeros(1,sizeT);
    idxR=zeros(1,sizeT);
    for t=1:1:sizeT
        ts=squeeze(profile(:,t));
        ts=interp1(1:1:length(ts),ts,1:interpStep:length(ts));
        idxsFrgrd=zeros(1,length(ts));
        [~,idxIni(t)]=max(ts);
        idxCur=idxIni(t);
        idxsFrgrd(idxCur)=1;
        
        idxL(t)=idxCur;
        idxR(t)=idxCur;
        
        stdFrgrd=std(ts(idxsFrgrd==1));
        stdBkgrd=std(ts(idxsFrgrd==0));
        
        stdSumIni=sqrt(sum(idxsFrgrd==1).*stdFrgrd.^2+sum(idxsFrgrd==0).*stdBkgrd.^2);
        stdSumCur=stdSumIni;
        stdSum2=stdSumCur;
        stdSumCur2=stdSum2;
        while stdSumCur2<=stdSum2 && idxL(t)>1 && idxR(t)<length(idxsFrgrd)
            stdSum2=stdSumCur2;
            stdSum=stdSumCur;
            idxCur=idxL(t);
            while stdSumCur<=stdSum && idxCur>1
                stdSum=stdSumCur;
                idxCur=idxCur-1;
                idxsFrgrd(idxCur)=1;
                stdFrgrd=std(ts(idxsFrgrd==1));
                stdBkgrd=std(ts(idxsFrgrd==0));
                stdSumCur=sqrt(sum(idxsFrgrd==1).*stdFrgrd.^2+sum(idxsFrgrd==0).*stdBkgrd.^2);
            end
            idxL(t)=idxCur;
            
            stdSum=stdSumCur;
            idxCur=idxR(t);
            while stdSumCur<=stdSum && idxCur<length(idxsFrgrd)
                stdSum=stdSumCur;
                idxCur=idxCur+1;
                idxsFrgrd(idxCur)=1;
                stdFrgrd=std(ts(idxsFrgrd==1));
                stdBkgrd=std(ts(idxsFrgrd==0));
                stdSumCur=sqrt(sum(idxsFrgrd==1).*stdFrgrd.^2+sum(idxsFrgrd==0).*stdBkgrd.^2);
            end
            idxR(t)=idxCur;
            stdSumCur2=stdSumCur;
        end
    end
    
    idxL=(idxL-1)*interpStep+1;
    idxR=(idxR-1)*interpStep+1;
    
end
%% Get the vessel bordere based on a drop to %.


%% Calculate the diameter and update the centerline
if medfiltS4>1
    idxL=medfilt1(idxL,medfiltS4,'omitnan','truncate');
    idxR=medfilt1(idxR,medfiltS4,'omitnan','truncate');
end
dP=idxR-idxL;

shiftC=round((round(idxL)+round(idxR))/2)-(floor(size(profile,1)/2)+1);

maskCenterLine=zeros(size(maskLine));
if (strcmp(chosenDirection,'H'))
    d=dP.*sind(abs(alpha));
    for t=1:1:sizeT
        indTS1=setTS1{t};
        indTS2=setTS2{t}+shiftC(t);
        indTS1(indTS2<1)=[];
        indTS2(indTS2<1)=[];
        indTS1(indTS2>sizeX)=[];
        indTS2(indTS2>sizeX)=[];
        
        idx1d=(indTS2-1).*sizeY+indTS1;
        tmpMask=zeros(sizeY,sizeX);
        tmpMask(idx1d)=1;
        maskCenterLine(:,:,t)=tmpMask;
        for i=1:1:floor((round(idxR)-round(idxL))/2)
            tmpMask=imdilate(tmpMask,SE);
        end
        maskVesselH(:,:,t)=tmpMask;
    end
else
    d=dP.*sind(90-abs(alpha));
    for t=1:1:sizeT
        indTS1=setTS1{t};
        indTS2=setTS2{t}+shiftC(t);
        indTS1(indTS2<1)=[];
        indTS2(indTS2<1)=[];
        indTS1(indTS2>sizeY)=[];
        indTS2(indTS2>sizeY)=[];
        
        idx1d=(indTS1-1).*sizeY+indTS2;
        tmpMask=zeros(sizeY,sizeX);
        tmpMask(idx1d)=1;
        maskCenterLine(:,:,t)=tmpMask;
        
        for i=1:1:floor((idxR-idxL)/2)
            tmpMask=imdilate(tmpMask,SE);
        end
        maskVesselV(:,:,t)=tmpMask;
    end
end


%% Assign results
if (strcmp(chosenDirection,'H'))
    region.maskVessel3=maskVesselH;
else
    region.maskVessel3=maskVesselV;
end

region.maskLine=maskLine;
region.profIdxL=idxL;
region.profIdxR=idxR;
region.maskCenterLine=maskCenterLine;
region.d=d;
region.dP=dP;
region.alpha=alpha;
region.profile=profile;

end