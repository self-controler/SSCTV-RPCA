function result=KRX(image,winlenth,targetShape)
%image 输入3维高光谱图像
%winlenth 局部窗口大小
%targetShape  目标形状
%result 检测结果
switch targetShape%目标形状
case 0
target_spatial=...
    [0 0 1 0 0;
    0 0 1 0 0;
    1 1 1 1 1;
    0 0 1 0 0;
    0 0 1 0 0];

case 1
target_spatial=...
    [1 1 1;
    1 1 1;
    1 1 1];
case 2 
target_spatial=...
    [1 0 0 0 1; 
    0 1 0 1 0;
    0 0 1 0 0;
    0 1 0 1 0;
    1 0 0 0 1];
case 3
target_spatial=...
    [1 0 0 0 0 0 1;
    0 1 0 0 0 1 0;
    0 0 1 0 1 0 0;
    0 0 0 1 0 0 0;
    0 0 1 0 1 0 0;
    0 1 0 0 0 1 0;
    1 0 0 0 0 0 1];
otherwise
disp('Unknown shape.')
end
sizeTarget=size(target_spatial);
% 若为二维图象，将其看作一幅第三维为1的三维图象
SizeA=size(image);
if length(SizeA)==2 
SizeA(3) = 1;
end
% 若为二维图象，将其看作一幅第三维为1的三维图象
SizeA=size(image);
if length(SizeA)==2 
SizeA(3) = 1;
end
%将image扩充为imageExpand
imageExpand((winlenth+1)/2:(winlenth-1)/2+SizeA(1), (winlenth+1)/2:(winlenth-1)/2+SizeA(2),:)=image; %中
imageExpand((winlenth-1)/2:-1:1,(winlenth+1)/2:(winlenth-1)/2+SizeA(2),:)=image(1:(winlenth-1)/2,:,:); %上
imageExpand(winlenth-1+SizeA(1):-1:(winlenth+1)/2+SizeA(1),(winlenth+1)/2:(winlenth-1)/2+SizeA(2),:)...%下
=image(SizeA(1)-(winlenth-3)/2:SizeA(1),:,:); %
imageExpand((winlenth+1)/2:(winlenth-1)/2+SizeA(1),(winlenth-1)/2:-1:1,:)=image(:,1:(winlenth-1)/2,:); %左
imageExpand((winlenth+1)/2:(winlenth-1)/2+SizeA(1),winlenth-1+SizeA(2):-1:(winlenth+1)/2+SizeA(2),:)...%右
=image(:,SizeA(2)-(winlenth-3)/2:SizeA(2),:); %
imageExpand((winlenth-1)/2:-1:1,(winlenth-1)/2:-1:1,:)=image(1:(winlenth-1)/2,1:(winlenth-1)/2,:); %左上
imageExpand(winlenth-1+SizeA(1):-1:(winlenth+1)/2+SizeA(1),winlenth-1+SizeA(2):-1:(winlenth+1)/2+... %右下
SizeA(2),:)=image(SizeA(1)-(winlenth-3)/2:SizeA(1),SizeA(2)-(winlenth-3)/2:SizeA(2),:); %
imageExpand(winlenth-1+SizeA(1):-1:(winlenth+1)/2+SizeA(1),(winlenth-1)/2:-1:1,:)... %左下
=image(SizeA(1)-(winlenth-3)/2:SizeA(1),1:(winlenth-1)/2,:); %
imageExpand((winlenth-1)/2:-1:1,winlenth-1+SizeA(2):-1:(winlenth+1)/2+SizeA(2),:)... %右上
=image(1:(winlenth-1)/2,SizeA(2)-(winlenth-3)/2:SizeA(2),:); %
clear image;

spatial_pattern=zeros(winlenth,winlenth);%winlenth应该是窗口大小
targetSize=size(target_spatial);
result=zeros(SizeA(1),SizeA(2));
for i=1:targetSize(1)
    for j=1:targetSize(2)
        spatial_pattern((winlenth-targetSize(1))/2+i,(winlenth-targetSize(1))/2+j)=1;
    end
end
s=zeros(1,winlenth*winlenth);
    k=1;
for i=1:winlenth
for j=1:winlenth

% s((i-1)*winlenth+j)=spatial_pattern(i,j);
s(k)=spatial_pattern(i,j);
k=k+1;
end
end
%%%%以上代码是将图像根据所需要检测的目标的形状进行扩充，变大
for i=(winlenth+1)/2:(winlenth-1)/2+SizeA(1)       %原图像区域
for j=(winlenth+1)/2:(winlenth-1)/2+SizeA(2)
dataWin=imageExpand(i-(winlenth-1)/2:i+(winlenth-1)/2,j-(winlenth-1)/2:j+(winlenth-1)/2,:);%从原始图像的第一个像素点开始遍历图像，截取winlength*winlength窗口大小立方体数据
x=[];
xkb=[];%存储外部窗口像素的光谱，用于计算Kb
xn=zeros(SizeA(3),1);
for m=1:winlenth
for n=1:winlenth %动
xn(:)=dataWin(m,n,:);     %取立方体的第一个像素点
x=[x,xn];   %%将每个像素点的数据排成一列，形成一个24（波段）*25的二维矩阵，行数为波段数，列数为像素数
if(m>3&&m<9&&n>3&&n<9);
else
    xkb=[xkb,xn];
end
end
end %窗 
%%%%%%以上是得到的样本矩阵代码，x样本数为winlenth*winlenth

%%%%%%%%计算KRX检测算子(Krt-Kut)'*Kb_inverse*(Krt-kut)
%步骤1：计算Kb_inverse
Kb=zeros(winlenth*winlenth-sizeTarget(1)*sizeTarget(2),winlenth*winlenth-sizeTarget(1)*sizeTarget(2));%未中心化的Kb
for k1=1:winlenth*winlenth-sizeTarget(1)*sizeTarget(2)
    for k2=1:k1
        Kb(k1,k2)=kernel(xkb(:,k1),xkb(:,k2));
    end
end
Kb=Kb+Kb';
for k=1:winlenth*winlenth-sizeTarget(1)*sizeTarget(2)
    Kb(k,k)=Kb(k,k)/2;
end
Im=(1.0/(winlenth*winlenth-sizeTarget(1)*sizeTarget(2)))*ones(winlenth*winlenth-sizeTarget(1)*sizeTarget(2),winlenth*winlenth-sizeTarget(1)*sizeTarget(2));
I=eye(winlenth*winlenth-sizeTarget(1)*sizeTarget(2));
Kb_center=(I-Im)*Kb*(I-Im);%中心化的Kb
beta=0.01;
%Kb_inverse=inv(Kb_center+beta*I);%计算Kb_center的逆Kb_inverse,正则化
Kb_inverse=pinv(Kb_center);
%[U,S,V]=svd(Kb_center);

%xx=x*s';
center=ceil(winlenth*winlenth/2);
xx=x(:,center);
%步骤2：计算Krt
kr=zeros(1,winlenth*winlenth-sizeTarget(1)*sizeTarget(2));
for k=1:winlenth*winlenth-sizeTarget(1)*sizeTarget(2)
    kr(1,k)=kernel(xx,xkb(:,k));
end
krmean=mean(kr);
Krt=kr-ones(size(kr))*krmean;

%步骤3：计算Kut
%ku=zeros(1,winlenth*winlenth-sizeTarget(1)*sizeTarget(2));
%kumean=0;
%for k1=1:winlenth*winlenth-sizeTarget(1)*sizeTarget(2)
   % for k2=1:winlenth*winlenth-sizeTarget(1)*sizeTarget(2)
      %  ku(1,k2)=ku(1,k2)+kernel(xkb(:,k1),xkb(:,k2));
    %end
%end
%for k=1:winlenth*winlenth-sizeTarget(1)*sizeTarget(2)
   % kumean=kumean+ku(1,k);
%end
ku=sum(Kb,1);
kumean=sum(ku);
kumean=kumean/((winlenth*winlenth-sizeTarget(1)*sizeTarget(2)).^2);
Kut=1.0/(winlenth*winlenth-sizeTarget(1)*sizeTarget(2))*ku-ones(size(ku))*kumean;

%步骤4：检测异常度
result(i-(winlenth-1)/2,j-(winlenth-1)/2)=   (Krt-Kut)*Kb_inverse*(Krt-Kut)';
end
end
BWImage=im2bw(result',0.42)%0.6524);

figure,imshow(result',[]);
figure,imshow(BWImage,[]);






























