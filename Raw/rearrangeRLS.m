%% rearrange RLS
fileName='D:\DLSI\data_4s_30mus_22800Hz.rls';
[data,sampling,timeStamps]=readRLS(fileName,0,0,[],'frame');

%modify data array
%data=data([1:800],:,:);
data=data(1:1024,:,:);

%save to RLS
fileWriteId = fopen(fileName,'w');

for k=1:(30*1024)
fwrite(fileWriteId,0 , 'uint64');
end
fseek(fileWriteId,0*1024,-1 );
fwrite(fileWriteId,uint64(size(data,1)), 'uint64');
fwrite(fileWriteId,uint64(size(data,2)), 'uint64');
fwrite(fileWriteId,uint64(size(data,3)), 'uint64');
fwrite(fileWriteId,uint64(sampling), 'uint64');

fseek(fileWriteId,30*1024,-1 );
for i=1:1:size(data,3)
    fwrite(fileWriteId,timeStamps(i), 'uint64');
    img=squeeze(data(:,:,i));
    fwrite(fileWriteId,img, 'uint8');    
    i
end
fclose(fileWriteId);

%%
fileName='D:\DLSI\Whole skull\wholeSkullDLSI.rls';
fileName2='E:\wholeSkullDLSIFix.rls';

fileReadId = fopen(fileName, 'r');
fileWriteId = fopen(fileName2,'w');

fseek(fileReadId,0*1024,-1 );
sizeX=fread(fileReadId,1,'*uint64');
sizeY=fread(fileReadId,1,'*uint64');
sizeT=fread(fileReadId,1,'*uint64');
sampling=fread(fileReadId,1,'*uint64');
fclose(fileReadId);

pts=[1:21];
for i=1:1:84
pts=[pts,i*30+1:i*30+1+21];
end
sizeX=length(pts);

for k=1:(30*1024)
fwrite(fileWriteId,0 , 'uint64');
end
fseek(fileWriteId,0*1024,-1 );
fwrite(fileWriteId,sizeX, 'uint64');
fwrite(fileWriteId,sizeY, 'uint64');
fwrite(fileWriteId,sizeT, 'uint64');
fwrite(fileWriteId,uint64(sampling), 'uint64');
fseek(fileWriteId,30*1024,-1 );

ii=0;
figure
while ii<sizeT
        ii
framesN=min(sizeT-ii,10000);
[data,sampling,timeStamps]=readRLS(fileName,ii,framesN,[],'frame');
ii=ii+framesN;

data=data-blackLevel;
data=data(pts,:,:);
imagesc(getSLSCI(squeeze(mean(data,3)),7,'gpu','none'));
pause(1)

for i=1:1:size(data,3)
    fwrite(fileWriteId,timeStamps(i), 'uint64');
    img=squeeze(data(:,:,i));
    fwrite(fileWriteId,img, 'uint8');   
end
end
fclose(fileWriteId);