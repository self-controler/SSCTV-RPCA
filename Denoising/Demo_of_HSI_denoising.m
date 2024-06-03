clear all;clc;
addpath(genpath('../SSCTV-RPCA-main/')) % windows platform
%% load data
hsi_name = 'pure_DCmall_small';
load([hsi_name,'.mat'])
clean_data=Ori_H;
clean_data       = Normalize(clean_data);
[M,N,p]        = size(clean_data);
tic
gaussian_level  = 0.0;
sparse_level   = 0.3;
noise_data       = GetNoise(clean_data,gaussian_level,sparse_level);
D = reshape(noise_data,[M*N,p]);
[mpsnr(1),mssim(1),ergas(1)]=msqia(clean_data, noise_data);

%% SSCTV-RPCA
it =2;
fprintf('======== SSCTV-RPCA  ========\n')
opts.rho = 1.1;
opts.lambda = 2/sqrt(M*N); 
tic;
[csstv_out,E] = ssctv_rpca(noise_data,opts);
t(it)=toc;
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, csstv_out);

%% CTV-RPCA
it =3;
fprintf('======== CTV-RPCA  ========\n')
opts.rho = 1.5;
opts.lambda = 3/sqrt(M*N); 
tic;
ctv_out = ctv_rpca(noise_data,opts);
t(it)=toc;
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, ctv_out);

%% RPCA
it =4;
D= zeros(M*N,p) ;
for ij=1:p
    bandp = noise_data(:,:,ij);
    D(:,ij)= bandp(:);
end
fprintf('========   RPCA  ========\n')
tic;
A_hat = rpca_m(D);
t(it)=toc;
rpca_out = reshape(A_hat,[M,N,p]);
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, rpca_out);

%% VBRPCA
it=5;
fprintf('========   VBRPCA  ========\n')
tic;
[X,~] = VBRPCA(D);
t(it)=toc;
vbrpca_out = reshape(X,[M,N,p]);
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, vbrpca_out);

%% LRTDTV
it=6;
fprintf('========   LRTDTV  ========\n')
tic;
[A_hat,~]=LRTDTV(noise_data,0.02,1,[100,100,20]);
t(it)=toc;
LRTDTV_out = reshape(A_hat,[M,N,p]);
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, LRTDTV_out);

%% LRTV
it =7;
fprintf('======== LRTV  ========\n')
tic;
[lrtv_out, out_value] = LRTV(noise_data,0.1,1,20);
t(it)=toc;
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, lrtv_out);

%% PRMF
it=8;
fprintf('========   PRMF  ========\n')
tic;
[U,V]=RPMF(D,5,1,1,1e-6);
t(it)=toc;
prmf_out = reshape(U*V,[M,N,p]);
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, prmf_out);

%% LRMR
it=9;
fprintf('========   LRMR  ========\n')
tic;
[A_hat,~]=MSILRMR(noise_data,0.01,10);
t(it)=toc;
LRMR_out = reshape(A_hat,[M,N,p]);
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, LRMR_out);

%% SSTV
it =10;
fprintf('======== SSTV  ========\n')
tic;
[sstv_out, out_value] = sstv(noise_data,0.3);
t(it)=toc;
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, sstv_out);

%% LRSSTV
it =11;
fprintf('======== LRSSTV  ========\n')
tic;
[lrsstv_out, out_value] = lrsstv(noise_data,0.2,2);
t(it)=toc;
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, lrsstv_out);
