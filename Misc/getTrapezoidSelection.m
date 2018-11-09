%getTrapezoidSelection - creates a trapezoid like selection data, where at
%least 2 lines are parallel to either X or Y axis.
%
% Syntax:  ROIs=getTrapezoidSelection(img)
%
% Inputs:
%    img  - image to select regions from
%
% Outputs:
%    ROIs - structure that contains trapezoid mask and bounding box
%           coordinates for the mask position in the original image
%
% Example:
%    ROIs=getTrapezoidSelection(img)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: DD Postnov, PhD
% BOAS lab, Boston University
% BMI, Copenhagen University
% email address: dpostnov@sund.ku.dk / dpostnov@bu.edu
% Last revision: 1-November-2018

%------------- BEGIN CODE --------------

function ROIs=getTrapezoidSelection(img)
disp('ROIs selection');
ii=1;
ifContinue=1;
while ifContinue==1
    title('Click on target points ')
    [x,y] = ginput(4);
    
    hold on
    plot(x,y,'k');
    hold off
    
    [~,idxy]=min(y);
    if idxy>1
        x=[x(idxy:end);x(1:idxy-1)];
        y=[y(idxy:end);y(1:idxy-1)];
    end
    l=zeros(1,length(x));
    for i=1:1:length(x)
        k=i+1;
        if k>length(x)
            k=1;
        end
        l(i)=sqrt((x(i)-x(k)).^2+(y(i)-y(k)).^2);
    end
    if l(1)+l(3)>l(2)+l(4)
        x=[x(1);x(end:-1:2)];
        y=[y(1);y(end:-1:2)];
    end
    
    m14=(y(1)-y(4))/(x(1)-x(4));
    b14=y(1)-m14*x(1);
    m23=(y(2)-y(3))/(x(2)-x(3));
    b23=y(2)-m23*x(2);
    
    hold on
    plot(x,y,'g');
    plot(x(1),y(1),'xg');
    plot(x(2),y(2),'og');
    hold off
    
    xx=x;
    yy=y;
    
    pointIdx=[1,1,2,2,3,3,4,4];
    
    xx(1)=x(2);
    yy(1)=m14*x(2)+b14;
    xx(2)=(y(2)-b14)./m14;
    yy(2)=y(2);
    xx(3)=x(1);
    yy(3)=m23*x(1)+b23;
    xx(4)=(y(1)-b23)./m23;
    yy(4)=y(1);
    xx(5)=x(4);
    yy(5)=m23*x(4)+b23;
    xx(6)=(y(4)-b23)./m23;
    yy(6)=y(4);
    xx(7)=x(3);
    yy(7)=m14*x(3)+b14;
    xx(8)=(y(3)-b14)./m14;
    yy(8)=y(3);
    
    dist=zeros(1,8);
    for i=1:1:8
        dist(i)=sqrt((xx(i)-x(pointIdx(i))).^2+(yy(i)-y(pointIdx(i))).^2);
    end
    
    [~, idx]=min(dist(pointIdx<3));
    
    
    x(pointIdx(idx))=xx(idx);
    y(pointIdx(idx))=yy(idx);
    
    if idx==1 || idx==3
        if dist(5)<dist(7)
            if xx(5)>1 && xx(5)<size(img,2) && yy(5)>1 && yy(5)<size(img,1)
                x(3)=xx(5);
                y(3)=yy(5);
            end
        else
            if xx(7)>1 && xx(7)<size(img,2) && yy(7)>1 && yy(7)<size(img,1)
                x(4)=xx(7);
                y(4)=yy(7);
            end
        end
    else
        if dist(6)<dist(8)
            if xx(6)>1 && xx(6)<size(img,2) && yy(6)>1 && yy(6)<size(img,1)
                x(3)=xx(6);
                y(3)=yy(6);
            end
        else
            if xx(8)>1 && xx(8)<size(img,2) && yy(8)>1 && yy(8)<size(img,1)
                x(4)=xx(8);
                y(4)=yy(8);
            end
        end
    end
    
    x(end+1)=x(1);
    y(end+1)=y(1);
    x=round(x);
    y=round(y);
    
    
    hold on
    plot(x,y,'r');
    hold off
    coordinates=cat(1,x,y);
    mask=poly2mask(coordinates(1:end/2),coordinates(end/2+1:end),size(img,1),size(img,2));
    maskPolyProps=regionprops(mask,'BoundingBox');
    maskCoord=maskPolyProps.BoundingBox;
    ROIs(ii).y=[floor(maskCoord(1))+1,floor(maskCoord(1)+maskCoord(3)-1)];
    ROIs(ii).x=[floor(maskCoord(2))+1,floor(maskCoord(2)+maskCoord(4)-1)];
    ROIs(ii).mask=mask(ROIs(ii).x(1):ROIs(ii).x(2),ROIs(ii).y(1):ROIs(ii).y(2));
    ii=ii+1;
    ifContinue = input('Select one more region? (1/0) ');
end
end