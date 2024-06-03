clc,clear
close all
addpath(genpath('ABU'));
addpath(genpath('proposed method'));
addpath(genpath('compare re'));
addpath(genpath('proposed result'));
addpath(genpath('comparison_methods'));
addpath(genpath('DATASET_NEW'));
k=1;

load RE1_NEW.mat
%% plane anomaly
if k==1;
    load Sandiego.mat
    M=Sandiego;   
    clear Sandiego;
    M(:,:,[1:6 33:35 97 104:110 153:166 221:224 ])=[];
    M(:,:,[94 95 96])=[];
    h=200;
    w=100;
    M=M(1:150,1:150,:);
    [no_lines,no_rows, no_bands]=size(M);
    sa=100;
    sb=100;
    load PlaneGT.mat 
    %0.9986  0.9549   0.9645(nm=2)
    GT=PlaneGT;

elseif k==2;
    load Urban.mat
    M=Urban; clear Urban
    M=M(1:80,189:288,:);
    [no_lines,no_rows, no_bands]=size(M);
    load UGt.mat
    GT=UGt;
    sa=80;
    sb=100;
%     w=80;  0.8910  0.9205(nm=2)
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

elseif k==4
        load abu-urban-2.mat
        M=data;
        [no_lines,no_rows, no_bands]=size(M);
        GT=map;
        sa=no_lines;
        sb=no_rows;
elseif k==5
        load abu-urban-1.mat
        M=data;
        [no_lines,no_rows, no_bands]=size(M);
        GT=map;
        sa=no_lines;
        sb=no_rows;
end

GT2=zeros(no_lines, no_rows);
if k~=2
GT2(1:sa,1:sb)=GT;
elseif k==2
    GT2(1:sa,end-sb+1:end)=GT;
end
M = M./max(M(:));
ndim = size(M);
M = hyperConvert2D(M);

GT=GT2;
%load re1.mat

% load d5_RE5
% re=RE5;
% thres = 0.05;
% set(gca,'Fontsize',20)

% [tpr,fpr,~] = roc(GT(:)',RE{1}(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
%plot(fpr_rx,tpr_rx,'-b','LineWidth',2)
re1 = reshape(RE{1},[no_lines,no_rows]);
figure
imagesc(re1)
axis image
axis off
axis normal;
set(gca,'position',[0 0 1 1]);
% [tpr,fpr,~] = roc(GT(:)',RE{2}(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'-k','LineWidth',1.5)
re2= reshape(RE{2},[no_lines,no_rows]);
figure
imagesc(re2)
axis image
axis off
axis normal;
set(gca,'position',[0 0 1 1]);

re3= reshape(RE{3},[no_lines,no_rows]);
figure
imagesc(re3)
axis image
axis off
axis normal;
set(gca,'position',[0 0 1 1]);

re4= reshape(RE{4},[no_lines,no_rows]);
figure
imagesc(re4)
axis image
axis off
axis normal;
set(gca,'position',[0 0 1 1]);

re5= reshape(RE{5},[no_lines,no_rows]);
figure
imagesc(re5)
axis image
axis off
axis normal;
set(gca,'position',[0 0 1 1]);

re6= reshape(RE{6},[no_lines,no_rows]);
figure
imagesc(re6)
axis image
axis off
axis normal;
set(gca,'position',[0 0 1 1]);

re7= reshape(RE{7},[no_lines,no_rows]);
figure
imagesc(re7)
axis image
axis off
axis normal;
set(gca,'position',[0 0 1 1]);

load D1_re_CP+TV+S_0.9967.mat
%re8 =reshape(sqrt(sum(E.^2)),no_lines,no_rows);
figure
imagesc(re)
axis image
axis off
axis normal;
set(gca,'position',[0 0 1 1]);

% [tpr,fpr,~] = roc(GT(:)',RE{3}(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'-y','LineWidth',1.5)
% 
% [tpr,fpr,~] = roc(GT(:)',RE{4}(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'-m','LineWidth',1.5)
% 
% [tpr,fpr,~] = roc(GT(:)',RE{5}(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'-c','LineWidth',1.5)
% 
% [tpr,fpr,~] = roc(GT(:)',RE{6}(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'-g','LineWidth',1.5)
% 
% [tpr,fpr,~] = roc(GT(:)',RE{7}(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'-b','LineWidth',1.5)
% 
% load E1_CPTVSNMF_0.9972.mat
% re8 =reshape(sqrt(sum(E.^2)),no_lines,no_rows);
% RE{8}=re8;
% [tpr,fpr,~] = roc(GT(:)',RE{8}(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'-r','LineWidth',1.5)
% legend('RX','RPCA','LRASR','GTVLRR','TRPCA','Degradation model 1','Degradation model 2','ATLRDSSU');
% xlabel('False alarm rate')
% ylabel('Probability of detection')
% set(gca,'XGrid','on');
% load D1_re1.mat
% [tpr,fpr,~] = roc(GT(:)',re1(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'Color',[0.5,0.16,0.16],'LineWidth',1.5)

% [tpr3,fpr3,thresholds3] = roc(GT(:)',re1(:)');
% AUC1=trapz(fpr3,tpr3);%figure
% 
% AUC=trapz(fpr,tpr);%figure
% hold on
% load D1_re2.mat
% [tpr,fpr,~] = roc(GT(:)',re2(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'-k','LineWidth',1.5)
% [tpr3,fpr3,thresholds3] = roc(GT(:)',re2(:)');
% AUC2=trapz(fpr3,tpr3);%figure
% 
% 
% load D1_re3.mat
% [tpr,fpr,~] = roc(GT(:)',re3(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'-y','LineWidth',1.5)
% 
% 
% [tpr3,fpr3,thresholds3] = roc(GT(:)',re3(:)');
% AUC3=trapz(fpr3,tpr3);%figure
% 
% load D1_re4.mat
% [tpr,fpr,~] = roc(GT(:)',re4(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'-m','LineWidth',1.5)
% 
% [tpr3,fpr3,thresholds3] = roc(GT(:)',re4(:)');
% AUC4=trapz(fpr3,tpr3);%figure
% 
% 
% load D1_re5.mat
% [tpr,fpr,~] = roc(GT(:)',re5(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'-c','LineWidth',1.5)
% [tpr3,fpr3,thresholds3] = roc(GT(:)',re5(:)');
% AUC5=trapz(fpr3,tpr3);%figure
% 
% 
% load D1_re6.mat
% [tpr,fpr,~] = roc(GT(:)',re6(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'-g','LineWidth',1.5)
% 
% [tpr3,fpr3,thresholds3] = roc(GT(:)',re6(:)');
% AUC6=trapz(fpr3,tpr3);%figure
% 
% load D1_re7.mat
% [tpr,fpr,~] = roc(GT(:)',re7(:)');
% tpr_rx = tpr(find(fpr<=thres));
% fpr_rx = fpr(find(fpr<=thres));
% plot(fpr_rx,tpr_rx,'-b','LineWidth',1.5)
% [tpr3,fpr3,thresholds3] = roc(GT(:)',re7(:)');
% AUC7=trapz(fpr3,tpr3);%figure
% 
% % load D1_re8.mat
% % [tpr,fpr,~] = roc(GT(:)',re8(:)');
% % tpr_rx = tpr(find(fpr<=thres));
% % fpr_rx = fpr(find(fpr<=thres));
% % %plot(fpr_rx,tpr_rx,'Color',[0.5,0.16,0.16],'LineWidth',1.5)
% % plot(fpr_rx,tpr_rx,'-r','LineWidth',1.5)
% 
% legend('RX','RPCA','LRASR','GTVLRR','TRPCA','SNMF','TV-SNMF','TV-CPSNMF');
% xlabel('False alarm rate')
% ylabel('Probability of detection')
% set(gca,'XGrid','on');
% figure
% plot(fpr1,tpr1,'Color',[0.5,0.16,0.16], 'LineWidth', 1.5);hold on;
% plot(fpr2,tpr2,'k-', 'LineWidth', 1.5);hold on;
% plot(fpr3,tpr3,'y-', 'LineWidth', 1.5);hold on;
% plot(fpr4,tpr4,'m-', 'LineWidth', 1.5);hold on;
% plot(fpr5,tpr5,'c-', 'LineWidth', 1.5);hold on;
% plot(fpr6,tpr6,'g-', 'LineWidth', 1.5);hold on;
% plot(fpr7,tpr7,'b-', 'LineWidth', 1.5);hold on;
% plot(fpr8,tpr8,'r-', 'LineWidth', 1.5);hold on;
% xlabel('false alarm rate')
% ylabel('probability of detection')
% legend('RX','RPCA','LRASR','GTVLRR','TRPCA','SNMF','TV-SNMF','TV-CPSNMF');
% hold off;
% %% plot ROC curve 
% figure
% plot(fpr1,tpr1,'b-', 'LineWidth', 1.5);hold on;
% plot(fpr2,tpr2,'k-', 'LineWidth', 1.5);hold on;
% plot(fpr3,tpr3,'r-', 'LineWidth', 1.5);hold on;
% plot(fpr4,tpr4,'y-', 'LineWidth', 1.5);hold on;
% plot(fpr5,tpr5,'m-', 'LineWidth', 1.5);hold on;
% plot(fpr6,tpr6,'k-', 'LineWidth', 1.5);hold on;
% plot(fpr7,tpr7,'c-', 'LineWidth', 1.5);hold on;
% plot(fpr8,tpr8,'g-', 'LineWidth', 1.5);hold on;
% 
% xlabel('false alarm rate')
% ylabel('probability of detection')
% legend('RPCA-RX','LRASR','RX','GTVLRR','TRPCA','SNMF','HSPRD','CP_TV_SNMF');
% hold off;














