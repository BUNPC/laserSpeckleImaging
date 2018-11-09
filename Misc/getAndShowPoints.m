function xy=getAndShowPoints(halfsize,colors,xy)
if isempty(xy)
count=1;
answ=1;
xy=[];

while answ==1
    [xx,yy]=ginput(1);
    hold on
    plot(xx,yy,'o','MarkerFaceColor',colors(count,:),'MarkerEdgeColor',[1,1,1],'MarkerSize',1.5*(halfsize*2+1));
    hold off
    xy(count,1)=int32(yy);
    xy(count,2)=int32(xx);
    count=count+1;
    answ=input('Select one more point (1/0)? ');
end
else
    for i=1:1:size(xy,1)
    hold on
    plot(xy(i,2),xy(i,1),'o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',[1,1,1],'MarkerSize',1.5*(halfsize*2+1));
    hold off   
    end
end
end