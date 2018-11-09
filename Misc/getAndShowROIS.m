function ROIs=getAndShowROIS(halfsize,colors,N,ROIs)
if isempty(ROIs)
count=1;
answ=1;
ROIs=[];

while answ==1
    
    BW = roipoly;
    s=regionprops(BW,'Orientation','PixelList');
    if length(s.PixelList)>N
    if abs(s.Orientation)<45
        s.PixelList=sortrows(s.PixelList,1);
        for i=N+1:1:length(s.PixelList)
            BW(s.PixelList(i,2),s.PixelList(i,1))=0;
        end
    else
        s.PixelList=sortrows(s.PixelList,2);
        for i=N+1:1:length(s.PixelList)
            BW(s.PixelList(i,2),s.PixelList(i,1))=0;
        end
        
    end
    
    
    hold on
    visboundaries(BW,'Color',colors(count,:));
    hold off
    ROIs(:,:,count)=BW;
    count=count+1;
    answ=input('Select one more point (1/0)? ');
    end
end
else
    for i=1:1:size(ROIs,3)
    BW = squeeze(ROIs(:,:,i));
    hold on
    visboundaries(BW,'Color',colors(i,:));
    hold off
    end
end
end