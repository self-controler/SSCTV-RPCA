clear all;clc;
% addpath(genpath('../CTV_code/')) % linux/MacOS platform
addpath(genpath('../SSCTV-main/')) % windows platform
for i=1:20
%% load data
hsi_name = 'pure_DCmall';
load([hsi_name,'.mat'])
if i==1
      clean_data       =  LoadCAVE('balloons_ms','data/CAVE/');
elseif i==2
      clean_data       =  LoadCAVE('beads_ms','data/CAVE/');  
elseif i==3
      clean_data       =  LoadCAVE('cd_ms','data/CAVE/'); 
elseif i==4
      clean_data       =  LoadCAVE('chart_and_stuffed_toy_ms','data/CAVE/'); 
elseif i==5
      clean_data       =  LoadCAVE('clay_ms','data/CAVE/'); 
elseif i==6
      clean_data       =  LoadCAVE('cloth_ms','data/CAVE/'); 
elseif i==7
      clean_data       =  LoadCAVE('egyptian_statue_ms','data/CAVE/');  
elseif i==8
      clean_data       =  LoadCAVE('face_ms','data/CAVE/'); 
elseif i==9
      clean_data       =  LoadCAVE('feathers_ms','data/CAVE/'); 
elseif i==10
      clean_data       =  LoadCAVE('flowers_ms','data/CAVE/'); 
elseif i==11
      clean_data       =  LoadCAVE('glass_tiles_ms','data/CAVE/');  
elseif i==12
      clean_data       =  LoadCAVE('hairs_ms','data/CAVE/'); 
elseif i==13
      clean_data       =  LoadCAVE('jelly_beans_ms','data/CAVE/');  
elseif i==14
      clean_data       =  LoadCAVE('oil_painting_ms','data/CAVE/');
elseif i==15
      clean_data       =  LoadCAVE('paints_ms','data/CAVE/');  
      elseif i==16
      clean_data       =  LoadCAVE('photo_and_face_ms','data/CAVE/');  
      elseif i==17
      clean_data       =  LoadCAVE('pompoms_ms','data/CAVE/');  
      elseif i==18
      clean_data       =  LoadCAVE('sponges_ms','data/CAVE/');  
      elseif i==19
      clean_data       =  LoadCAVE('superballs_ms','data/CAVE/');  
      elseif i==20
      clean_data       =  LoadCAVE('thread_spools_ms','data/CAVE/');  
end
clean_data=Ori_H;
clean_data       = Normalize(clean_data);
[M,N,p]        = size(clean_data);
tic
gaussian_level = 0;
sparse_level   = 0.3;
noise_data       = GetNoise(clean_data,gaussian_level,sparse_level);
D = reshape(noise_data,[M*N,p]);
mpsnr = zeros(3,1);
mssim = zeros(3,1);
ergas = zeros(3,1);
[mpsnr(1),mssim(1),ergas(1)]=msqia(clean_data, noise_data);


%% SSCTV-RPCA
it =2;
fprintf('======== CSSTV-RPCA  ========\n')
opts.rho = 1.1;
opts.lambda = 2/sqrt(M*N); 
tic;
[csstv_out,E] = ssctv_rpca(noise_data,opts);
t(it)=toc;
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, csstv_out);
%% CTV-RPCA
it =3;
fprintf('======== CTV-RPCA  ========\n')
opts.rho = 1.25;
opts.lambda =3/sqrt(M*N); 
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
fprintf('======== LRTV  ========\n')
tic;
[sstv_out, out_value] = sstv(noise_data,0.3);
t(it)=toc;
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, sstv_out);
%% LRSSTV
it =11;
fprintf('======== LRTV  ========\n')
tic;
[lrsstv_out, out_value] = lrsstv(noise_data,0.2,2);
t(it)=toc;
[mpsnr(it),mssim(it),ergas(it)]=msqia(clean_data, lrsstv_out);
ps(i,:)=mpsnr;
ss(i,:)=mssim;
er(i,:)=ergas;
end

%showband = 13;
%figure;
%Y = clean_data(:,:,showband);
%subplot(2,2,1);imshow(Y,[]);title('original band')
%Y = rpca_out(:,:,showband);
%subplot(2,2,2);imshow(Y,[]);title(['rpca, psnr:',num2str(mpsnr(4))])
%Y = ctv_out(:,:,showband);
%subplot(2,2,3);imshow(Y,[]);title(['ctv-rpca, psnr:',num2str(mpsnr(2))])
%Y = csstv_out(:,:,showband);
%subplot(2,2,4);imshow(Y,[]);title(['csstv-rpca, psnr:',num2str(mpsnr(3))])