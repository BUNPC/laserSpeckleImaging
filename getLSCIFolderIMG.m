%getLSCIFolderIMG - Reads raw burst files recorded by speckleSoftware 
%(A. Dunn) in the specified folder and processes them in chosen type of 
%laser speckle contrast
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3,input4,input5,input6,input7)
%
% Inputs:
%    folderName - folder name as string
%    burstStart - first burst to read
%    burstN     - number of bursts to read, or [] if read all bursts
%                 starting from burstStart
%    lscType    - type of laser speckle contraast: use 'sLSCI' or 'tLSCI'
%    kernelSize - number of pixels in the kernel (or in the side of it)
%    procType   - choose the processor type: use 'cpu' or 'gpu'
%    dsType     - downsampling type result is either same size as data or 
%                 downsampled by kernel size. Use: 'none' or 'kernel'
%
% Outputs:
%    LSCI         - processed laser speckle data as [y,x,t] 3d matrix
%    time         - time in seconds, 0 for the first burst.
%
% Example:
%    [LSCI,time]=getLSCIFolderIMG('D:\trial 1 left',1,100,'sLSCI',7,'gpu','full');
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

function [LSCI,time]=getLSCIFolderIMG(folderName,burstStart,burstN,lscType,kernelSize,procType,dsType)
workingDir=pwd;
cd(folderName);
dataFiles=dir('img.*');
dataFiles(count({dataFiles.name},'.')>1)=[]; % exclude all files except img stacks
[~,order]=sort({dataFiles.name});
dataFiles=dataFiles(order);

if ~isempty(burstN)
    dataFiles=dataFiles(burstStart:burstStart+burstN-1);
end
[subdata, ~, timeStamp]=readSingleIMG(dataFiles(1).name);

imgsPerBurst=size(subdata,3);
filesN=length(dataFiles);

if strcmp(lscType,'tLSCI') && imgsPerBurst<kernelSize
    ME = MException('Error:WrongParams', ...
        'Number of frames in a burst is less than tLSCI kernelSize');
    throw(ME);
end

time=zeros(1,filesN);
time(1)=getTime(timeStamp);

if strcmp(lscType,'sLSCI')
    subLSCI=getSLSCI(subdata,kernelSize,procType,dsType);
elseif strcmp(lscType,'tLSCI')
    subLSCI=getTLSCI(subdata,kernelSize,procType,'kernel');
end
subLSCI=mean(subLSCI,3);

LSCI=zeros(size(subLSCI,1),size(subLSCI,2),filesN,'single');
LSCI(:,:,1)=subLSCI;
for i=2:1:filesN
    [subdata, ~, timeStamp]=readSingleIMG(dataFiles(i).name);
    
    if strcmp(lscType,'sLSCI')
        subLSCI=getSLSCI(subdata,kernelSize,procType,dsType);
    elseif strcmp(lscType,'tLSCI')
        subLSCI=getTLSCI(subdata,kernelSize,procType,'kernel');
    end
    subLSCI=mean(subLSCI,3);
    LSCI(:,:,i)=subLSCI;
    time(i)=getTime(timeStamp);
end
time=time-time(1);

cd(workingDir);

    function t=getTime(tstamp)
        %Does not work if recording spans 2 different months.
        t=str2double(tstamp(15:17))./1000+str2double(tstamp(13:14))+str2double(tstamp(11:12))*60+...
            str2double(tstamp(9:10))*60*60+str2double(tstamp(7:8))*60*60*24;
    end

end

%------------- END OF CODE --------------
% Comments: Speckle contrast calculated from images in the same burst is
% averaged. Time calculation will fail in case if recording spans 2
% different months.
