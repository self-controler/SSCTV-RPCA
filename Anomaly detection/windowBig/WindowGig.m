function Y = WindowGig(X,local,Osz,times,outX)
%  ������ͼ��ȡһ�����ڷŴ�ĳ���
%  X������ͼ����[0,1]֮���ֵ����ɫ��ڰ׵ģ�Y�������ͼ��һ���ǲ�ɫ�ģ���һ����ά��������
%  local ��Ҫ�Ŵ�Ŀ�����Ͻǵ�λ��,��ռԭͼ�ı�������ʾ
%  Osz   �ǿ��ԭ��С����ռԭͼ�ı�������ʾ
%  times ��Ҫ�Ŵ�ı���
%  outX  �ǷŴ�Ŀ����ͼ��Ķ�С����0.1�����������ԭͼ��С��0.1��

if nargin<2
    local  = [0.2,0.21];
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
    error('��������ȷ��ͼ���ʽ')
end
Patch   = Y(x1(1):x2(1),x1(2):x2(2),1:3);
rePatch = imresize(Patch,times);

resizeP = size(rePatch);
y2      = sizeY(1:2);
y1      = max(1,y2-resizeP(1:2)+1);
Y(y1(1):y2(1),y1(2):y2(2),1:3)=rePatch;

%ԭ��߿���
for i = 1:3
    if i==2
        Y(x1(1):x2(1),[x1(2),x2(2)],i)=1;
        Y([x1(1),x2(1)],x1(2):x2(2),i)=1;
    else
        Y(x1(1):x2(1),[x1(2),x2(2)],i)=0;
        Y([x1(1),x2(1)],x1(2):x2(2),i)=0;
    end
    if i==1
        Y(y1(1):y2(1),[y1(2),y2(2)],i)=1;
        Y([y1(1),y2(1)],y1(2):y2(2),i)=1;
    else
        Y(y1(1):y2(1),[y1(2),y2(2)],i)=0;
        Y([y1(1),y2(1)],y1(2):y2(2),i)=0;
    end    
end

for i = [2,3]
    
end
end