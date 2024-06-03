function result=KRX(image,winlenth,targetShape)
%image ����3ά�߹���ͼ��
%winlenth �ֲ����ڴ�С
%targetShape  Ŀ����״
%result �����
switch targetShape%Ŀ����״
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
% ��Ϊ��άͼ�󣬽��俴��һ������άΪ1����άͼ��
SizeA=size(image);
if length(SizeA)==2 
SizeA(3) = 1;
end
% ��Ϊ��άͼ�󣬽��俴��һ������άΪ1����άͼ��
SizeA=size(image);
if length(SizeA)==2 
SizeA(3) = 1;
end
%��image����ΪimageExpand
imageExpand((winlenth+1)/2:(winlenth-1)/2+SizeA(1), (winlenth+1)/2:(winlenth-1)/2+SizeA(2),:)=image; %��
imageExpand((winlenth-1)/2:-1:1,(winlenth+1)/2:(winlenth-1)/2+SizeA(2),:)=image(1:(winlenth-1)/2,:,:); %��
imageExpand(winlenth-1+SizeA(1):-1:(winlenth+1)/2+SizeA(1),(winlenth+1)/2:(winlenth-1)/2+SizeA(2),:)...%��
=image(SizeA(1)-(winlenth-3)/2:SizeA(1),:,:); %
imageExpand((winlenth+1)/2:(winlenth-1)/2+SizeA(1),(winlenth-1)/2:-1:1,:)=image(:,1:(winlenth-1)/2,:); %��
imageExpand((winlenth+1)/2:(winlenth-1)/2+SizeA(1),winlenth-1+SizeA(2):-1:(winlenth+1)/2+SizeA(2),:)...%��
=image(:,SizeA(2)-(winlenth-3)/2:SizeA(2),:); %
imageExpand((winlenth-1)/2:-1:1,(winlenth-1)/2:-1:1,:)=image(1:(winlenth-1)/2,1:(winlenth-1)/2,:); %����
imageExpand(winlenth-1+SizeA(1):-1:(winlenth+1)/2+SizeA(1),winlenth-1+SizeA(2):-1:(winlenth+1)/2+... %����
SizeA(2),:)=image(SizeA(1)-(winlenth-3)/2:SizeA(1),SizeA(2)-(winlenth-3)/2:SizeA(2),:); %
imageExpand(winlenth-1+SizeA(1):-1:(winlenth+1)/2+SizeA(1),(winlenth-1)/2:-1:1,:)... %����
=image(SizeA(1)-(winlenth-3)/2:SizeA(1),1:(winlenth-1)/2,:); %
imageExpand((winlenth-1)/2:-1:1,winlenth-1+SizeA(2):-1:(winlenth+1)/2+SizeA(2),:)... %����
=image(1:(winlenth-1)/2,SizeA(2)-(winlenth-3)/2:SizeA(2),:); %
clear image;

spatial_pattern=zeros(winlenth,winlenth);%winlenthӦ���Ǵ��ڴ�С
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
%%%%���ϴ����ǽ�ͼ���������Ҫ����Ŀ�����״�������䣬���
for i=(winlenth+1)/2:(winlenth-1)/2+SizeA(1)       %ԭͼ������
for j=(winlenth+1)/2:(winlenth-1)/2+SizeA(2)
dataWin=imageExpand(i-(winlenth-1)/2:i+(winlenth-1)/2,j-(winlenth-1)/2:j+(winlenth-1)/2,:);%��ԭʼͼ��ĵ�һ�����ص㿪ʼ����ͼ�񣬽�ȡwinlength*winlength���ڴ�С����������
x=[];
xkb=[];%�洢�ⲿ�������صĹ��ף����ڼ���Kb
xn=zeros(SizeA(3),1);
for m=1:winlenth
for n=1:winlenth %��
xn(:)=dataWin(m,n,:);     %ȡ������ĵ�һ�����ص�
x=[x,xn];   %%��ÿ�����ص�������ų�һ�У��γ�һ��24�����Σ�*25�Ķ�ά��������Ϊ������������Ϊ������
if(m>3&&m<9&&n>3&&n<9);
else
    xkb=[xkb,xn];
end
end
end %�� 
%%%%%%�����ǵõ�������������룬x������Ϊwinlenth*winlenth

%%%%%%%%����KRX�������(Krt-Kut)'*Kb_inverse*(Krt-kut)
%����1������Kb_inverse
Kb=zeros(winlenth*winlenth-sizeTarget(1)*sizeTarget(2),winlenth*winlenth-sizeTarget(1)*sizeTarget(2));%δ���Ļ���Kb
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
Kb_center=(I-Im)*Kb*(I-Im);%���Ļ���Kb
beta=0.01;
%Kb_inverse=inv(Kb_center+beta*I);%����Kb_center����Kb_inverse,����
Kb_inverse=pinv(Kb_center);
%[U,S,V]=svd(Kb_center);

%xx=x*s';
center=ceil(winlenth*winlenth/2);
xx=x(:,center);
%����2������Krt
kr=zeros(1,winlenth*winlenth-sizeTarget(1)*sizeTarget(2));
for k=1:winlenth*winlenth-sizeTarget(1)*sizeTarget(2)
    kr(1,k)=kernel(xx,xkb(:,k));
end
krmean=mean(kr);
Krt=kr-ones(size(kr))*krmean;

%����3������Kut
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

%����4������쳣��
result(i-(winlenth-1)/2,j-(winlenth-1)/2)=   (Krt-Kut)*Kb_inverse*(Krt-Kut)';
end
end
BWImage=im2bw(result',0.42)%0.6524);

figure,imshow(result',[]);
figure,imshow(BWImage,[]);






























