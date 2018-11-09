%readFolderIMG - Reads all raw files recorded by speckleSoftware (A. Dunn)
%found in the folder
%
% Syntax:  [output1,output2,output3,output4] = function_name(input1,input2,input3)
%
% Inputs:
%    folderName   - folder name as string
%    burstStart   - first burst to read
%    burstN       - number of bursts to read, or [] if read all bursts
%                   starting from burstStart
%
% Outputs:
%    data         - raw speckle data as [y,x,t] 3d matrix
%    imgsPerBurst - number of images per burst
%    exposure     - exposure time in microseconds
%    timeStamp    - times when bursts were generated as yyyymmddxxxxxxxx
%
% Example:
%    [data,imgsPerBurst,exposure,timeStamps]=readFolderIMG('D:\Occlusion',1,[])
%
% Other m-files required: readSingleIMG.m
% Subfunctions: none
% MAT-files required: none
%
% See also: readSingleIMG.m

% Author: DD Postnov, PhD
% BOAS lab, Boston University
% email address: dpostnov@bu.edu
% Last revision: 3-March-2018

%------------- BEGIN CODE --------------

function [data,imgsPerBurst,exposure,timeStamps]=readFolderIMG(folderName,burstStart,burstN)
workingDir=pwd;
cd(folderName);
dataFiles=dir('img.*');
dataFiles(count({dataFiles.name},'.')>1)=[]; % exclude all files except img stacks
[~,order]=sort({dataFiles.name});
dataFiles=dataFiles(order);

if ~isempty(burstN)
dataFiles=dataFiles(burstStart:burstStart+burstN-1);
end

[subdata, exposure, timeStamp]=readSingleIMG(dataFiles(1).name);

sizeY=size(subdata,1);
sizeX=size(subdata,2);
imgsPerBurst=size(subdata,3);
filesN=length(dataFiles);
data=zeros(sizeY,sizeX,imgsPerBurst.*filesN,'uint8');
timeStamps=cell(filesN,1);

data(:,:,1:imgsPerBurst)=subdata;
timeStamps{1}=timeStamp;

for i=2:1:filesN
    [subdata, exposure, timeStamp]=readSingleIMG(dataFiles(i).name);
    data(:,:,(i-1)*imgsPerBurst+1:i*imgsPerBurst)=uint8(subdata);
    timeStamps{i}=timeStamp;
end

cd(workingDir);
end

%------------- END OF CODE --------------
% Comments: recording protocol of the speckleSoftware assumes burst
% capture mode, which leads to a varying sampling time between frames in
% output data matrix. Also it is assumed that raw data was saved in uint8
% format.
