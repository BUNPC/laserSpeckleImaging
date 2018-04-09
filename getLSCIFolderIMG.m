%getLSCIFolderIMG - Reads raw burst files recorded by speckleSoftware
%(A. Dunn) in the specified folder and processes them in chosen type of
%laser speckle contrast
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3,input4,input5,input6,input7)
%
% Inputs:
%    procType     - run as single thread or as multicore: use 'cluster' or
%                   'thread'. Cluster is default, compatible only with CPU
%                   laser speckle processing
%    folderName   - folder name as string
%    burstStart   - first burst to read
%    burstN       - number of bursts to read, or [] to read all bursts
%                   starting from burstStart
%    lscType      - type of laser speckle contraast: use 'sLSCI' or 'tLSCI'
%                   sLSCI is default.
%    kernelSize   - number of pixels in the kernel (or in the side of it)
%    procTypeLSCI - choose the processor type: use 'cpu' or 'gpu'.
%                   'cpu' is default
%    dsType       - downsampling type result is either same size as data or
%                   downsampled by kernel size. Use: 'none' or 'kernel'.
%                   'none' is default.
%    toCrop       - Crop flag. If 1, then allows you to select the crop
%                   area before processing.
%
% Outputs:
%    LSCI         - processed laser speckle data as [y,x,t] 3d matrix
%    time         - time in seconds, first value counted from the first day
%                   of the month.
%
% Example:
%    [LSCI,time]=getLSCIFolderIMG('cluster','D:\trial 1 left',1,[],'sLSCI',7,'cpu','full');
%
% Other m-files required: readSingleIMG.m, getSLSCI.m, getTLSCI.m
% Subfunctions: none
% MAT-files required: none
%
% See also: readSingleIMG.m, readFolderIMG.m, getSLSCI.m, getTLSCI.m

% Author: DD Postnov, PhD
% BOAS lab, Boston University
% email address: dpostnov@bu.edu
% Last revision: 3-March-2018

%------------- BEGIN CODE --------------

function [LSCI,time,imgsPerBurst]=getLSCIFolderIMG(procType,folderName,burstStart,burstN,lscType,kernelSize,procTypeLSCI,dsType,toCrop)
% add directory with code to path
[workDir]=fileparts(mfilename('fullpath'));
addpath(workDir);

if strcmp(procTypeLSCI,'cluster') && strcmp(procType,'cluster')
    ME = MException('Error:WrongParams', ...
        'Execution of cluster inside of cluster');
    throw(ME);
end

% find raw files and select the subset of interest
cd(folderName);
dataFiles=dir('img.*');
dataFiles(count({dataFiles.name},'.')>1)=[];
[~,order]=sort({dataFiles.name});
dataFiles=dataFiles(order);
if ~isempty(burstN)
    dataFiles=dataFiles(burstStart:burstStart+burstN-1);
end

% get and process the first file to fetch data dimensions and time estimate
tic
[subdata, ~, timeStamp]=readSingleIMG(dataFiles(1).name);
elapsedTime=toc;

if toCrop==1
    figure
    imagesc(squeeze(mean(subdata,3)));
    title('Select crop area')
    h=imrect;
    pos = wait(h);
    pos=round(pos);
    cropY=pos(2):1:pos(2)+pos(4);
    cropX=pos(1):1:pos(1)+pos(3);
    subdata=subdata(cropY,cropX,:);
end

imgsPerBurst=size(subdata,3);
filesN=length(dataFiles);

if strcmp(lscType,'tLSCI') && imgsPerBurst<kernelSize
    ME = MException('Error:WrongParams', ...
        'Number of frames in a burst is less than tLSCI kernelSize');
    throw(ME);
end

if strcmp(lscType,'tLSCI')
    subLSCI=getTLSCI(subdata,kernelSize,procTypeLSCI,dsType);
else %lscType='sLSCI'
    subLSCI=getSLSCI(subdata,kernelSize,procTypeLSCI,dsType);
end
subLSCI=mean(subLSCI,3);
elapsedTime=elapsedTime+toc;

LSCI=zeros(size(subLSCI,1),size(subLSCI,2),filesN,'single');
time=zeros(1,filesN);

% process the data
disp(['Current time is ', datestr(now,'HH:mm:ss')]);
if strcmp(procType,'thread')
    time(1)=getTime(timeStamp);
    LSCI(:,:,1)=subLSCI;
    disp([num2str(filesN),' burst files to process. Time remaining ~',num2str(elapsedTime*(filesN-1)),'s'])
    for i=2:1:filesN
        [subdata, ~, timeStamp]=readSingleIMG(dataFiles(i).name);
        if toCrop==1
            subdata=subdata(cropY,cropX,:);
        end
        if strcmp(lscType,'TLSCI')
            subLSCI=getTLSCI(subdata,kernelSize,procTypeLSCI,dsType);
        else %lscType='sLSCI'
            subLSCI=getSLSCI(subdata,kernelSize,procTypeLSCI,dsType);
        end
        subLSCI=mean(subLSCI,3);
        LSCI(:,:,i)=subLSCI;
        time(i)=getTime(timeStamp);
    end
else %cluster case
    pool=gcp;
    tic
    parfor i=1:1:pool.NumWorkers
        [subdata, ~, timeStamp]=readSingleIMG(dataFiles(i).name);
        if toCrop==1
            subdata=subdata(cropY,cropX,:);
        end
        if strcmp(lscType,'tLSCI')
            subLSCI=getTLSCI(subdata,kernelSize,procTypeLSCI,dsType);
        else %lscType='sLSCI'
            subLSCI=getSLSCI(subdata,kernelSize,procTypeLSCI,dsType);
        end
        subLSCI=mean(subLSCI,3);
        LSCI(:,:,i)=subLSCI;
        timeStamps{i}=timeStamp;
    end
    elapsedTime=toc;
    elapsedTime=elapsedTime./pool.NumWorkers;
    disp([num2str(filesN),' burst files to process. Time remaining ~',num2str(elapsedTime*(filesN-pool.NumWorkers)),'s'])
    
    parfor i=pool.NumWorkers+1:1:filesN
        [subdata, ~, timeStamp]=readSingleIMG(dataFiles(i).name);
        if toCrop==1
            subdata=subdata(cropY,cropX,:);
        end
        if strcmp(lscType,'tLSCI')
            subLSCI=getTLSCI(subdata,kernelSize,procTypeLSCI,dsType);
        else %lscType='sLSCI'
            subLSCI=getSLSCI(subdata,kernelSize,procTypeLSCI,dsType);
        end
        subLSCI=mean(subLSCI,3);
        LSCI(:,:,i)=subLSCI;
        timeStamps{i}=timeStamp;
    end
    for i=1:1:filesN
        time(i)=getTime(timeStamps{i});
    end
end

cd(workDir);

    function t=getTime(tstamp)
        %Does not work if recording spans 2 different months.
        t=str2double(tstamp(15:17))./1000+str2double(tstamp(13:14))+str2double(tstamp(11:12))*60+...
            str2double(tstamp(9:10))*60*60+str2double(tstamp(7:8))*60*60*24;
    end

end

%------------- END OF CODE --------------
% Comments: Speckle contrast calculated from images in the same burst is
% averaged. Time values calculation will fail in case if recording spans 2
% different months.
