clear all;clc
%addpath(genpath('../../CTV_code/')) % linux/MacOS platform
addpath(genpath('..\..\CTV_code\')) % windows platform
h = 20;
w = 20;
band = 200;
r_s  = 0.2;
rho_s  = 0.2;
smooth_flag = 1;
smo=1;
for rhot=1:60
    rho_s=0.01*rhot;
[D,RLmat,RSmat]=generate_M(h,w,band,r_s,rho_s,smooth_flag);
diff_v = diff3(D,[h,w,band]);
normD = norm(D,'fro');
noise_data = reshape(D,[h,w,band]);
%% CSSTV-RPCA
out = csstv_rpca(noise_data);
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
out = csstv_rpca(noise_data,opts);
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
cssr(rhot)=log10(err_csstv);
cr(rhot)=log10(err_ctv);
rpr(rhot)=log10(err_rpca);
end
a=0.01:0.01:0.6;
plot(a,cssr,'ro-',a,cr,'b.-',a,rpr,'go-');
xlabel('\rho_s');
ylabel('log_{10}||X_0-X||/||X_0||');
legend('csstv','ctv','rpca','Location','Best');
%fprintf('======== Result ============\n')
%fprintf(' error of ssctv-rpca :  %.6f\n ',err_csstv);
%fprintf('error of   ctv-rpca :  %.6f\n ',err_ctv);
%fprintf('error of       rpca :  %.6f\n ',err_rpca);

