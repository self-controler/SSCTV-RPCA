clear all;clc;
% addpath(genpath('../CTV_code/')) % linux/MacOS platform
addpath(genpath('../ccstv_junheng/')) % windows platform
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


