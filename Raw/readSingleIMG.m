%readSingleIMG - Reads one burst/raw file recorded by speckleSoftware (A. Dunn)
%
% Syntax:  [output1,output2,output3] = function_name(input1)
%
% Inputs:
%    fileName  - file name as string
%
% Outputs:
%    data      - raw speckle data as [y,x,t] 3d matrix
%    exposure  - exposure time in microseconds
%    timeStamp - time when file was generated as yyyymmddhhmmssxxx, 
%                xxx are ms
%
% Example: 
%    [data, exposure, timeStamp]=readSingleIMG('D:\img.20171205145924109')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: readFolderIMG.m

% Author: DD Postnov, PhD
% BOAS lab, Boston University
% email address: dpostnov@bu.edu
% Last revision: 3-March-2018

%------------- BEGIN CODE --------------

function [data,exposure,timeStamp]=readSingleIMG(fileName)
  timeStamp=strsplit(fileName,'.');
  timeStamp=timeStamp{2};
  fileReadId=fopen(fileName,'r');
  sizeX=fread(fileReadId,1,'unsigned short');
  sizeY=fread(fileReadId,1,'unsigned short');
  sizeT=fread(fileReadId,1,'unsigned short');
  exposure=uint64(fread(fileReadId,1,'unsigned short')*1000); %convert to microseconds
  data=fread(fileReadId,sizeY*sizeX*sizeT,'uchar');
  data=reshape(data,[sizeY sizeX sizeT]);
  fclose(fileReadId);
end

%------------- END OF CODE --------------
% Comments:
