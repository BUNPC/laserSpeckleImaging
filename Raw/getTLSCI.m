%getTLSCI - calculates temporal Laser Speckle Contrast Images
%
% Syntax:  output1 = function_name(input1,input2,input3,input4)
%
% Inputs:
%    data       - raw laser speckle data as 3d [y,x,t] matrix
%    kernelSize - number of pixels in a side of the kernel
%    procType   - choose the processor type: use 'cpu' or 'gpu'
%                 'cpu' is default
%    dsType     - downsampling type result is either same size as data or
%                 downsampled by kernel size. Use: 'none' or 'kernel'
%                 'none' is default
%
% Outputs:
%    tLSCI      - processed data as [y,x,t] 3d matrix
%
% Example:
%    tLSCI=getSLSCI(data,25,'gpu','none')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getSLSCI.m

% Author: DD Postnov, PhD
% BOAS lab, Boston University
% email address: dpostnov@bu.edu
% Last revision: 3-March-2018

%------------- BEGIN CODE --------------

function tLSCI=getTLSCI(data,kernelSize,procType,dsType)
if strcmp(dsType,'kernel')
    T=1:kernelSize:size(data,3)-kernelSize+1;
else  % dsType='none'
    T=1:1:size(data,3)-kernelSize+1;
end
tLSCI=zeros(size(data,1),size(data,2),length(T),'single');
counter=1;
if strcmp(procType,'gpu')
    for i=T
        frames=gpuArray(single(data(:,:,i:i+kernelSize-1)));
        frameMean=squeeze(mean(frames,3));
        frameSTD=squeeze(std(frames,0,3));
        tLSCI(:,:,counter)=gather(frameSTD./frameMean);
        counter=counter+1;
    end
else %procType='cpu' - single thread cpu processing
    for i=T
        frames=single(data(:,:,i:i+kernelSize-1));
        frameMean=squeeze(mean(frames,3));
        frameSTD=squeeze(std(frames,0,3));
        tLSCI(:,:,counter)=frameSTD./frameMean;
        counter=counter+1;
    end
end
end

%------------- END OF CODE --------------
% Comments: large input data can lead to the memory overflow, particularily
% when uint8 imput data is provided. This can be controlled by additional
% outer loop and/or by conversion sLSCI data to scaled integer or by
% allowing downsampling by kernel size



