clear;
clc

addpath(genpath('comparison_methods'));
addpath(genpath('dataset'));
addpath(genpath('prox_operator'));
addpath(genpath('tensor_toolbox'));
addpath(genpath('tensor_SVD'));
addpath(genpath('utils'));

Method_list = {'RX', 'RPCA-RX', 'LRASR', 'LSMAD', 'GTVLRR', 'TRPCA', 'PTA','RPCA','CTV','CSSTV'};
% run or not
Run_RX      = 1;
Run_RPCA_RX = 1;
Run_LRASR   = 1;
Run_LSMAD   = 1;
Run_GTVLRR  = 1;
Run_TRPCA   = 1;
Run_PTA     = 1;
Run_RPCA   = 1;
Run_ctv   = 1;
Run_ssctv   = 1;
en_list = [];
%% plane anomaly
k=2;
if k==1
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
elseif k==5
    load Salians_syn.mat
    M=hsi;
    clear hsi;
    [no_lines, no_rows, no_bands]=size(M);
    load Salians_syn_gt.mat
    GT = hsi_gt;
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
X=X';
X=reshape(X,no_lines,no_rows,no_bands);

%% show data and groundtruth of anomaly map
figure
image = zeros(no_lines,no_rows,3);
image(:,:,1) = M_3D(:,:,49);
image(:,:,2) = M_3D(:,:,27);
image(:,:,3) = M_3D(:,:,7);
maxP = max(image(:));minP = min(image(:));
image = (image-minP)/(maxP-minP);


%% RX 
if Run_RX == 1
    disp(['Running', Method_list{1},'...']);
    re{1} = reshape(RxDetector1(M),[no_lines,no_rows]);
    [tpr{1},fpr{1},~] = roc(GT2(:)',re{1}(:)');
    AUC(1)=trapz(fpr{1},tpr{1});
    fprintf('AUC of RX:%f\n',AUC(1));
    en_list = [en_list, 1];
end

%% RPCA-RX
if Run_RPCA_RX == 1
    disp(['Running', Method_list{2},'...']);
    lambda = 0.006;
    [r0 ,Output_S, ~] = Unsupervised_RPCA_Detect_v1(X,lambda);
    XS = reshape(Output_S, no_lines*no_rows, no_bands);
    r = RX(XS'); 
    re{2} = reshape(r,[no_lines,no_rows]);
    [tpr{2},fpr{2},~] = roc(GT2(:)',re{2}(:)');
    AUC(2)=trapz(fpr{2},tpr{2});
    fprintf('AUC of RPCA_RX:%f\n',AUC(2));
    en_list = [en_list, 2];
end

%% LRASR
if Run_LRASR == 1
    beta=0.1;lamda=0.1;
    Y = M;
    Dict=ConstructionD_lilu(Y,15,20);
    [~,E]=LRASR(Y,Dict,beta,lamda,1);
    re{3} =reshape(sqrt(sum(E.^2)),no_lines,no_rows);
    [tpr{3},fpr{3},~] = roc(GT2(:)',re{3}(:)');
    AUC(3)=trapz(fpr{3},tpr{3});
    fprintf('AUC of LRASR:%f\n',AUC(3));
    en_list = [en_list, 3];
end

%% LSMAD
if Run_LSMAD == 1
    Y = M;
    [L,~,RMSE,error]=GoDec(Y',28,floor(0.0022*no_bands)*9,2);
    L=L';
    mu=mean(L,2);
    r_new2=(diag((Y-repmat(mu,[1,no_lines*no_rows]))'*pinv(cov(L'))*(Y-repmat(mu,[1,no_lines*no_rows]))))';
    re{4} = reshape(r_new2,[no_lines,no_rows]);
    [tpr{4},fpr{4},~] = roc(GT2(:)',re{4}(:)');
    AUC(4)=trapz(fpr{4},tpr{4});
    fprintf('AUC of LSMAD:%f\n',AUC(4));
    en_list = [en_list, 4];
end

%% GTVLRR
if Run_GTVLRR == 1
    lambda = 0.5;
    beta = 0.5;
    gamma =0.01;
    display = true;
    [~,S] = lrr_tv_manifold(Y,Dict,lambda,beta,gamma,[no_lines,no_rows],display);
    S = reshape(S',[no_lines,no_rows,no_bands]);
    re{5}=reshape(sqrt(sum(S.^2,3)),[no_lines,no_rows]);
    [tpr{5},fpr{5},thresholds5] = roc(GT2(:)',re{5}(:)');
    AUC(5)=trapz(fpr{5},tpr{5});
    fprintf('AUC of GTVLRR:%f\n',AUC(5));
    en_list = [en_list, 5];
end

%% TRPCA
if Run_TRPCA == 1
    %lambda = 0.01;
    lambda = 1/sqrt(200*186);
    [~,E] = trpca_tnn(M_3D,lambda);
    E = reshape(E,no_lines,no_rows,no_bands);
    re{6}=reshape(sqrt(sum(E.^2,3)),[no_lines,no_rows]);
    [tpr{6},fpr{6},thresholds6] = roc(GT2(:)',re{6}(:)');
    AUC(6)=trapz(fpr{6},tpr{6});%figure
    fprintf('AUC of TRPCA:%f\n',AUC(6));
    en_list = [en_list, 6];
end

%% PTA
if Run_PTA == 1
    addpath(genpath('C:\Users\LENOVO\Desktop\办公文件\Anomaly detection'));
    %%%%
    mask_reshape = reshape(GT2, 1, no_lines*no_rows);
    anomaly_map = logical(double(mask_reshape)>0);
    normal_map = logical(double(mask_reshape)==0);
    %%%%
    %%%%PTA with LTV-norm
    tol1=1e-4;
    tol2=1e-6;
    maxiter=400;
    truncate_rank=1;
    alphia=1.7;
    beta=0.069;
    tau=0.1;
    
    tic;
    [X,S,area] = AD_Tensor_LILU1(M_3D,alphia,beta,tau,truncate_rank,maxiter,tol1,tol2,normal_map,anomaly_map);
    toc
    re{7}=reshape(sqrt(sum(S.^2,3)),[no_lines,no_rows]);
    [tpr{7},fpr{7},~] = roc(GT2(:)',re{7}(:)');
    AUC(7)=trapz(fpr{7},tpr{7});%figure
    fprintf('AUC of PTA:%f\n',AUC(7));
    en_list = [en_list, 7];
end
%% RPCA
%lambda = 0.01;
if Run_RPCA == 1
    lambda = 1/sqrt(no_rows*no_bands);
    [~,E,~] = rpca_m(M');
    E = reshape(E,no_lines,no_rows,no_bands);
    re{8}=reshape(sqrt(sum(E.^2,3)),[no_lines,no_rows]);
    [tpr{8},fpr{8},thresholds8] = roc(GT2(:)',re{8}(:)');
    AUC(8)=trapz(fpr{8},tpr{8});%figure
    fprintf('AUC of TRPCA:%f\n',AUC(8));
    en_list = [en_list,8];
end 
%% CTV
%lambda = 0.01;
if Run_ctv == 1
    opts.lambda = 1/sqrt(no_rows*no_bands);
    [~,E] = ctv_rpca(M_3D,opts);
    E = reshape(E,no_lines,no_rows,no_bands);
    re{9}=reshape(sqrt(sum(E.^2,3)),[no_lines,no_rows]);
    [tpr{9},fpr{9},thresholds9] = roc(GT2(:)',re{9}(:)');
    AUC(9)=trapz(fpr{9},tpr{9});%figure
    fprintf('AUC of TRPCA:%f\n',AUC(9));
    en_list = [en_list,9];
end
%% SSCTV
%lambda = 0.01;
if Run_ssctv == 1
opts.lambda = 1/sqrt(no_rows*no_bands);
    [~,E] = ssctv_rpca(M_3D,opts);
    E = reshape(E,no_lines,no_rows,no_bands);
    re{10}=reshape(sqrt(sum(E.^2,3)),[no_lines,no_rows]);
    [tpr{10},fpr{10},thresholds10] = roc(GT2(:)',re{10}(:)');
    AUC(10)=trapz(fpr{10},tpr{10});%figure
    fprintf('AUC of TRPCA:%f\n',AUC(10));
    en_list = [en_list,10];
end





%% plot ROC curve 
figure
for i = 1:length(Method_list)
    subplot(2,5,i); imagesc(re{i});axis off;
end



