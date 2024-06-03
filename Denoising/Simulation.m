clear all;clc
%addpath(genpath('../../CTV_code/')) % linux/MacOS platform
addpath(genpath('..\..\Denoising\')) % windows platform
h = 20;
w = 20;
band = 200;
r_s  = 0.2;%default
rho_s  = 0.2;%default
smooth_flag = 1;
cssm=zeros(50,50);
cm=zeros(50,50);
rpcam=zeros(50,50);
for rm=1:50
    for rhom=1:50
        r_s  = 0.01*rm;
        rho_s  =0.01*rhom ;   
        sumcss=0;
        sumc=0;
        sumr=0;
        for t=1:20

[D,RLmat,RSmat]=generate_M(h,w,band,r_s,rho_s,smooth_flag);
diff_v = diff3(D,[h,w,band]);
normD = norm(D,'fro');
noise_data = reshape(D,[h,w,band]);
%% SSCTV-RPCA
out = ssctv_rpca(noise_data);
A_hat = reshape(out,[h*w,band]);
it = 1;
err(it) = norm(A_hat-RLmat,'fro')/normD;
%% CTV-RPCA
out = ctv_rpca(noise_data);
A_hat = reshape(out,[h*w,band]);
it = 2;
err(it) = norm(A_hat-RLmat,'fro')/normD;
%% RPCA
A_hat = rpca_m(D);
it = 3;
err(it) = norm(A_hat-RLmat,'fro')/normD;
%% CSSTV with varing lambda
ratio = max(2-rho_s*2,1+r_s*2); % need to fine tune for CSSTV
clear opts;
opts.lambda = 2*ratio/sqrt(h*w);
out = ssctv_rpca(noise_data,opts);
A_hat = reshape(out,[h*w,band]);
it = 4;
err(it) = norm(A_hat-RLmat,'fro')/normD;
%% CTV with varing lambda
ratio = max(2-rho_s*2,1+r_s*2);
clear opts;
opts.lambda = 3*ratio/sqrt(h*w);
out = ctv_rpca(noise_data,opts);
A_hat = reshape(out,[h*w,band]);
it = 5;
err(it) = norm(A_hat-RLmat,'fro')/normD;
%% RPCA with varing lambda
lambda = ratio/sqrt(h*w);
A_hat = rpca_m(D,lambda);
it = 6;
err(it) = norm(A_hat-RLmat,'fro')/normD;
err_csstv = min(err(1),err(4));
err_ctv = min(err(2),err(5));
err_rpca = min(err(3),err(6));
if err_csstv<0.01
    sumcss=sumcss+1;
end
if err_ctv<0.01
    sumc=sumc+1;
end
if err_rpca<0.01
    sumr=sumr+1;
end
        end
        cssm(rm,rhom)=sumcss/20;
        cm(rm,rhom)=sumc/20;
        rpcam(rm,rhom)=sumr/20;
        if sumcss==0&&sumc==0&&sumr==0
            break
        end
    end
end
%fprintf('======== Result ============\n')
%fprintf(' error of ssctv-rpca :  %.6f\n ',err_csstv);
%fprintf('error of   ctv-rpca :  %.6f\n ',err_ctv);
%fprintf('error of       rpca :  %.6f\n ',err_rpca);