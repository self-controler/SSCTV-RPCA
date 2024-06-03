clear;
clc

addpath(genpath('comparison_methods'));
addpath(genpath('dataset'));
addpath(genpath('tensorlab_2016-03-28'));
addpath(genpath('prox_operator'));
addpath(genpath('tensor_toolbox'));
addpath(genpath('tensor_SVD'));
k=1;
%% plane anomaly
if k==1;
    load Sandiego.mat
    M=Sandiego;
    
    clear Sandiego;
    M(:,:,[1:6 33:35 97 104:110 153:166 221:224 ])=[];
    M(:,:,[94 95 96])=[];
    h=200;
    w=100;
    % M=M(265:364,1:100,:);
    M=M(1:150,1:150,:);
    [no_lines,no_rows, no_bands]=size(M);
    sa=100;
    sb=100;
    
    load PlaneGT.mat 
    %0.9986  0.9549
    GT=PlaneGT;

elseif k==2;
    load Urban.mat
    M=Urban; clear Urban
    M=M(1:80,189:288,:);
%     M=M(1:150,139:288,:);
    [no_lines,no_rows, no_bands]=size(M);
    load UGt.mat
    GT=UGt;
%     GT(70:end, 1:10)=0;
%     GT(70:end,90:end)=0;
    sa=80;
    sb=100;
%     h=100;
%     w=80;  0.8910
elseif k==3;
    load Sandiego.mat
    M=Sandiego;
    
    clear Sandiego;
    M(:,:,[1:6 33:35 97 104:110 153:166 221:224 ])=[];
    M(:,:,[94 95 96])=[];
    M=M(111:210,181:380,:);
     [no_lines,no_rows, no_bands]=size(M);
    sa=100;
    sb=200;
    load PlaneGT2.mat

%     GT=PlaneGT;0.9488  % 0.9212

elseif k==4
        load RIT_F1.mat
        M=PM;
        [no_lines,no_rows, no_bands]=size(M);
        GT=GTf1;
        clear Mf1 GTf1
        DicT=sig1;
        [a b]=meshgrid(40:20:160);
        tau=[0.2 0.3 0.4  0.5 0.6 0.7 0.8];
        for i=1:7
            temp=zeros(1,5,no_bands);
            for j=1:7
                temp(:,j,:)=DicT;
            end
            M(a(i,i),b(:,i),:)=(1-tau(i))*M(a(i,i),b(:,i),:)+tau(i)*temp;
        end
        GT=zeros(size(GT));
        for i=1:7
            
            for j=1:7
                GT(a(i,j),b(i,j))=1;
            end
            
        end%%end simulated
        
        sa=no_lines;
        sb=no_rows;
end

%M_3D= M./max(M(:));
ndim = size(M);
M = hyperConvert2D(M);
M=M./repmat(sqrt(sum(M.^2)),size(M,1),1);
M_3D= hyperConvert3D(M,no_lines,no_rows);
GT2=zeros(no_lines, no_rows);
X = M;
if k~=2
GT2(1:sa,1:sb)=GT;
elseif k==2
    GT2(1:sa,end-sb+1:end)=GT;
end

% K = 15;%cluster center
% P = 20;%number of pixels in each cluster
% % Dictionary construction
% Dic = DicCon_KSAC(X, K, P);
% solve representation coefficient matrix S and saprse error matrix E
X=X';
X=reshape(X,no_lines,no_rows,no_bands);

% %% RX 0.964452
% re1 = RxDetector1(M);
% [tpr1,fpr1,thresholds1] = roc(GT2(:)',re1(:)');
% AUC1=trapz(fpr1,tpr1);
% re1=reshape(sqrt(sum(re1.^2,3)),[no_lines,no_rows]);
% figure
% imagesc(re1)
% axis off
% fprintf('AUC RX:%f\n',AUC1);

%% RPCA-RX
lambda = 0.006;
[r0 ,Output_S, Output_L] = Unsupervised_RPCA_Detect_v1(X,lambda);
XS = reshape(Output_S, no_lines*no_rows, no_bands);
r = RX(XS'); 
re2 = reshape(r,[no_lines,no_rows]);
figure
imagesc(re2)
axis image
title('RPCA-RX for simulated data')
[tpr2,fpr2,thresholds2] = roc(GT2(:)',re2(:)');
AUC2=trapz(fpr2,tpr2);%figure
fprintf('AUC:%f\n',AUC2);

%% LRASR
beta=0.1;
lamda=0.1;
 Y = M;
 Dict=ConstructionD_lilu(Y,15,20);
[S,E]=LRASR(Y,Dict,beta,lamda,1);
re3 =reshape(sqrt(sum(E.^2)),no_lines,no_rows);
figure
imagesc(re3)
axis image
title('LRASR for simulated data')
[tpr3,fpr3,thresholds3] = roc(GT2(:)',re3(:)');
AUC3=trapz(fpr3,tpr3);%figure
fprintf('AUC:%f\n',AUC3);

%% LSMAD
 Y = M;
[L,S,RMSE,error]=GoDec(Y',28,floor(0.0022*no_bands)*9,2);
L=L';
S=S';
mu=mean(L,2);
r_new2=(diag((Y-repmat(mu,[1,no_lines*no_rows]))'*pinv(cov(L'))*(Y-repmat(mu,[1,no_lines*no_rows]))))';
re4 = reshape(r_new2,[no_lines,no_rows]);
figure
imagesc(re4)
axis image
title('LSMAD for simulated data')
[tpr4,fpr4,thresholds4] = roc(GT2(:)',re4(:)');
AUC4=trapz(fpr4,tpr4);%figure
fprintf('AUC:%f\n',AUC4);

%% GTVLRR
lambda = 0.5;
beta = 0.5;
gamma =0.001;
display = true;
[Q,S] = lrr_tv_manifold(Y,Dict,lambda,beta,gamma,[no_lines,no_rows],display);
re5=sqrt(sum(S.^2));
figure
imagesc(reshape(re5,[no_lines,no_rows]))
axis image
title('GTVLRR for simulated data')
[tpr5,fpr5,thresholds5] = roc(GT2(:)',re5(:)');
AUC5=trapz(fpr5,tpr5);%figure
fprintf('AUC:%f\n',AUC5);

%% TRPCA
lambda = 0.01;
% addpath(genpath('tensorlab_2016-03-28'));
% addpath(genpath('prox_operator'));
% addpath(genpath('tensor_toolbox'));
% addpath(genpath('tensor_SVD'));
[L,E] = trpca_tnn(M_3D,lambda);
E = reshape(E,no_lines,no_rows,no_bands);
re6=reshape(sqrt(sum(E.^2,3)),[no_lines,no_rows]);
figure
imagesc(re6)
axis image
title('TPCA for simulated data')
[tpr6,fpr6,thresholds6] = roc(GT2(:)',re6(:)');
AUC6=trapz(fpr6,tpr6);%figure
fprintf('AUC:%f\n',AUC6);

load Sandiego.mat
imshow(reshape(M(10,:),no_lines,no_rows),[]);
title('Band 10 image')
figure(2)
%% plot ROC curve 
figure
plot(fpr1,tpr1,'b-', 'LineWidth', 1.5);hold on;
plot(fpr2,tpr2,'k-', 'LineWidth', 1.5);hold on;
plot(fpr3,tpr3,'r-', 'LineWidth', 1.5);hold on;
plot(fpr4,tpr4,'y-', 'LineWidth', 1.5);hold on;
plot(fpr5,tpr5,'m-', 'LineWidth', 1.5);hold on;
plot(fpr6,tpr6,'k-', 'LineWidth', 1.5);hold on;
% plot(fpr7,tpr7,'c-', 'LineWidth', 1.5);hold on;
% plot(fpr8,tpr8,'g-', 'LineWidth', 1.5);hold on;

xlabel('false alarm rate')
ylabel('probability of detection')
legend('RX','RPCA-RX','LRASR','LSMAD','GTVLRR','TPCA');

hold off;
