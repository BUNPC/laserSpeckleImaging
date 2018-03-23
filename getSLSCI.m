%getSLSCI - calculates spatial Laser Speckle Contrast Images using square
%kernel
%
% Syntax:  output1 = function_name(input1,input2,input3,input4)
%
% Inputs:
%    data       - raw laser speckle data as 3d [y,x,t] matrix
%    kernelSize - number of pixels in a side of the kernel
%    procType   - choose the processor type: use 'cpu' or 'gpu'
%    dsType     - downsampling type result is either full frame or frame 
%                 downsampled by kernel size. Use: 'none' or 'kernel'
%
% Outputs:
%    sLSCI      - processed data as [y,x,t] 3d matrix
%
% Example: 
%    sLSCI=getSLSCI(data,7,'gpu','none')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getTLSCI.m

% Author: DD Postnov, PhD
% BOAS lab, Boston University
% email address: dpostnov@bu.edu
% Last revision: 3-March-2018

%------------- BEGIN CODE --------------

function sLSCI=getSLSCI(data,kernelSize,procType,dsType)
if strcmp(dsType,'none')
Y=1:1:size(data,1);
X=1:1:size(data,2);
elseif strcmp(dsType,'kernel')
halfSize=floor(kernelSize/2);
Y=halfSize+1:kernelSize:size(data,1)-halfSize;
X=halfSize+1:kernelSize:size(data,2)-halfSize;
else
Y=1:1:size(data,1);
X=1:1:size(data,2);
end

sLSCI=zeros(length(Y),length(X),size(data,3),'single');

if strcmp(procType,'cpu')
    for i=1:1:size(data,3)
        frame=single(data(:,:,i));
        frameMean=imfilter(frame,fspecial('average',[kernelSize kernelSize]));
        frameSTD=stdfilt(frame,ones(kernelSize));
        sLSCI(:,:,i)=frameSTD(Y,X)./frameMean(Y,X);
    end
elseif strcmp(procType,'gpu')
    for i=1:1:size(data,3)
        frame=gpuArray(single(data(:,:,i)));
        frameMean=imfilter(frame,fspecial('average',[kernelSize kernelSize]));
        frameSTD=stdfilt(frame,ones(kernelSize));
        sLSCI(:,:,i)=gather(frameSTD(Y,X)./frameMean(Y,X));
    end
end
end

%------------- END OF CODE --------------
% Comments: large input data can lead to the memory overflow, particularily
% when uint8 imput data is provided. This can be controlled by additional
% outer loop and/or by conversion sLSCI data to scaled integer or by
% allowing downsampling by kernel size