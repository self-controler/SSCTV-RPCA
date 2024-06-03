function Y = WindowGig_left(X,local,Osz,times,outX)
%  用来将图像取一个窗口放大的程序
%  X是输入图像是[0,1]之间的值，彩色或黑白的，Y，是输出图像，一定是彩色的，是一个三维的张量。
%  local 是要放大的块的左上角的位置,以占原图的比例来表示
%  Osz   是块的原大小，以占原图的比例来表示
%  times 是要放大的倍数
%  outX  是放大的块伸出图像的多小，如0.1代表伸出数是原图大小的0.1倍

if nargin<2
    local  = [0.1,0.11];
end
if nargin<3
    Osz    = [0.2,0.2];
end
if nargin<4
    times  = 2;
end
if nargin<5
    outX   = 0.1;
end

sizeX   = size(X);
sizeP   = round(sizeX(1:2).*Osz);
sizeOut = round(sizeX(1:2)*outX);
sizeY   = [sizeX(1:2)+sizeOut,3];

x1      = max(round(sizeX(1:2).*local),1);
x2      = min(x1+sizeP,sizeX(1:2));

Y       = ones(sizeY);
if length(sizeX)==2
    for i = 1:3
        Y(1:sizeX(1),1:sizeX(2),i) = X;
    end
elseif length(sizeX)==3
    Y(1:sizeX(1),1:sizeX(2),1:3) = X;
else
    error('请输入正确的图像格式')
end
Patch   = Y(x1(1):x2(1),x1(2):x2(2),1:3);
rePatch = imresize(Patch,times);

resizeP = size(rePatch);
y2      = sizeY(1:2);
y1      = max(1,y2-resizeP(1:2)+1);
Y(y1(1):y2(1),1:y2(2)-y1(2)+1,1:3)=rePatch;

%原块边框线
for i = 1:3
    if i==2
        Y(x1(1):x2(1),[x1(2),x2(2)],i)=1;
        Y([x1(1),x2(1)],x1(2):x2(2),i)=1;
    else
        Y(x1(1):x2(1),[x1(2),x2(2)],i)=0;
        Y([x1(1),x2(1)],x1(2):x2(2),i)=0;
    end
    if i==2 || i==1
        Y(y1(1):y2(1),[1,y2(2)-y1(2)+1],i)=1;
        Y([y1(1),y2(1)],1:y2(2)-y1(2)+1,i)=1;
    else
        Y(y1(1):y2(1),[1,y2(2)-y1(2)+1],i)=0;
        Y([y1(1),y2(1)],1:y2(2)-y1(2)+1,i)=0;
    end    
end

for i = [2,3]
    
end
end