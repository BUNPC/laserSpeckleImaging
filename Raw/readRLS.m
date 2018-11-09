function [data,sampling,timeStamps]=readRLS(fileName,startT,sizeT,ROI,type)
fileReadId = fopen(fileName, 'r');
fseek(fileReadId,0*1024,-1 );
%NOTE IN SOME RLS FILES Y and X were misplaced!!! Check the image and fix
%it!
sizeX=fread(fileReadId,1,'*uint64');
sizeY=fread(fileReadId,1,'*uint64');

if sizeT==0
    sizeT=fread(fileReadId,1,'*uint64');
else
    fread(fileReadId,1,'*uint64');
end
sampling=fread(fileReadId,1,'*uint64');
timeStamps=zeros(sizeT,1,'int64');
firstByte=30*1024+sizeX*sizeY*startT+8*startT;
fseek(fileReadId,firstByte,-1 );

if strcmp(type,'frame')
    data=zeros(sizeX,sizeY,sizeT,'uint8');
    for t=1:1:sizeT
        timeStamps(t)=fread(fileReadId,1,'*uint64');
        data(:,:,t)=fread(fileReadId,[sizeX,sizeY],'*uint8');
    end
elseif strcmp(type,'rect')
    data=zeros(length(ROI(1,1):1:ROI(1,2)),length(ROI(2,1):1:ROI(2,2)),sizeT,'uint8');
    for t=1:1:sizeT
        timeStamps=fread(fileReadId,1,'*uint64');
        frame=fread(fileReadId,[sizeX,sizeY],'*uint8');
        data(:,:,t)=frame(ROI(1,1):1:ROI(1,2),ROI(2,1):1:ROI(2,2));
    end
elseif strcmp(type,'1d')
    data=zeros(length(ROI),sizeT,'uint8');
    for t=1:1:sizeT
        timeStamps=fread(fileReadId,1,'*uint64');
        frame=fread(fileReadId,[sizeX,sizeY],'*uint8');
        data(:,t)=frame(ROI);
    end
end
fclose(fileReadId);

end