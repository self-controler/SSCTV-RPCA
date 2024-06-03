%% solve the problem 
%              min_X ||D_zD_x(X)||_*+||D_zD_y(X)||_*  +\lambda||E||_1
%                                   s.t.  Y= X+E
%                          ===============================
%              min_X ||X31||_*+||X32||_* +\lambda||E||_1
%                            s.t.  Y= X+E
%                                  D_x(X)=X1 
%                                  D_y(X)=X2 
%                                  D_z(X1)=X31
%                                  D_z(X2)=X32
%                          ===============================                       
%         D is difference operator,T is difference tensor,T is known
%  ------------------------------------------------------------------------


function [output_image,E] = ssctv_rpca(noise_data,opts)
[M,N,p] = size(noise_data);
if ~exist('opts','var')
    opts=[]; 
end
if isfield(opts,'maxIter')
    maxIter = opts.maxIter; 
else
    maxIter = 200; 
end
if isfield(opts,'rho')
    rho = opts.rho; 
else
    rho = 1.03; 
end
if isfield(opts,'tol')
    tol = opts.tol; 
else
    tol = 1e-6; 
end
if isfield(opts,'lambda')
    lambda = opts.lambda; 
else
    lambda = 2/sqrt(M*N); 
end
sizeD   = size(noise_data);
D       = zeros(M*N,p) ;
for i=1:p
    bandp = noise_data(:,:,i);
    D(:,i)= bandp(:);
end
normD   = norm(D,'fro');
% initialize
norm_two = lansvd(D, 1, 'L');
norm_inf = norm( D(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);

mu  = 1/dual_norm;%1.25/norm_two % this one can be tuned
mu1 = 1*mu;
max_mu = mu * 1e7;
%% FFT setting
h               = sizeD(1);
w               = sizeD(2);
d               = sizeD(3);
Eny_x   = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ;
Eny_y   = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;
determ  = Eny_x + Eny_y;

Eny_z   = ( abs(psf2otf([+1, -1], [w,d,h])) ).^2  ;
Eny_z   =  permute(Eny_z, [3, 1 2]);
determz =  Eny_z;
%% Initializing optimization variables
X              = D;
G1             = zeros(M*N,p);%reshape(diff_x(X,sizeD),sizeD);
G2             = zeros(M*N,p);%reshape(diff_y(X,sizeD),sizeD);
E              = zeros(M*N,p);
%M1 =zeros(size(D));  % multiplier for D-X-E
M1 = D / dual_norm;
M2 = M1;%zeros(size(D));  % multiplier for Dx_X-G1
M3 = M2;%zeros(size(D));  % multiplier for Dy_X-G2
M4 = M3;%zeros(size(D));  % multiplier for Dz_G1-G31
M5 = M4;%zeros(size(D));  % multiplier for Dz_G2-G32
% main loop
iter = 0;
tic
while iter<maxIter
    iter          = iter + 1;   
    %% -Updata G32,G31
    [u,s,v] = svd(reshape(diff_z(G1,sizeD),[M*N,p])+M4/mu1,'econ');
    G31     = u*softthre(s,1/mu1)*v';
    [u,s,v] = svd(reshape(diff_z(G2,sizeD),[M*N,p])+M5/mu1,'econ');
    G32     = u*softthre(s,1/mu1)*v';
    %% -Updata G1,G2
    diffT_p = diff_zT(mu1*G31 - M4,sizeD);
    numer1  = reshape(diffT_p + mu1*diff_x(X,sizeD) + M2(:),sizeD);
    x       = real( ifftn( fftn(numer1) ./ (mu1*determz + mu1) ) );
    G1      = reshape(x,[M*N,p]);
    diffT_p = diff_zT(mu1*G32 - M5,sizeD);
    numer1  = reshape(diffT_p + mu1*diff_y(X,sizeD) + M3(:),sizeD);
    x       = real( ifftn( fftn(numer1) ./ (mu1*determz + mu1) ) );
    G2      = reshape(x,[M*N,p]);
    
    %% -Updata X
    diffT_p  = diff_xT(mu1*G1-M2,sizeD)+diff_yT(mu1*G2-M3,sizeD);
    numer1   = reshape( diffT_p + mu*(D(:)-E(:)) + M1(:), sizeD);
    x        = real( ifftn( fftn(numer1) ./ (mu1*determ + mu ) ) );
    X        = reshape(x,[M*N,p]);
    %% -Update E
    E             = softthre(D-X+M1/mu, lambda/mu);
    %% stop criterion  
    leq1 = D -X -E;
    leq2 = reshape(diff_x(X,sizeD),[M*N,p])- G1;
    leq3 = reshape(diff_y(X,sizeD),[M*N,p])- G2;
    leq4 = reshape(diff_z(G1,sizeD),[M*N,p])- G31;
    leq5 = reshape(diff_z(G2,sizeD),[M*N,p])- G32;
    
    stopC1 = norm(leq1,'fro')/normD;
    stopC2 = max(abs(leq2(:)));
    stopC4 = norm(leq4,'fro')/normD;
    if mod(iter,10)==0
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e')  ...
                ',Y-X-E=' num2str(stopC1,'%2.3e') ',||DX-X1||=' num2str(stopC2,'%2.3e')...
                ',|DZ-X3|' num2str(stopC4,'%2.3e')]);
    end
    if stopC1<tol && stopC2<tol
        break;
    else
        M1  = M1 + mu*leq1;
        M2  = M2 + mu1*leq2;
        M3  = M3 + mu1*leq3;
        M4  = M4 + mu1*leq4;
        M5  = M5 + mu1*leq5;
        mu  = min(max_mu,mu*rho); 
        mu1 = min(max_mu,mu1*rho); 
    end 
end
output_image = reshape(X,[M,N,p]);
end