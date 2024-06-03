function [Xten,Sten] = lrsstv(Yten,lambda,beta)
sizeD =size(Yten);
mu = 0.01;
rho = 1.25;
maxIter = 200;
tol =1e-5;
normD   = norm(Yten(:),'fro');
%% FFT setting
h               = sizeD(1);
w               = sizeD(2);
d               = sizeD(3);
Eny_x   = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ;
Eny_y   = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;
determ  = Eny_x + Eny_y;

Eny_z   = ( abs(psf2otf([+1, -1], [w,d,h]))).^2  ;
Eny_z   =  permute(Eny_z, [3, 1 2]);

Lten = zeros(h,w,d);
Xten = zeros(h,w,d);
Sten = zeros(h,w,d);
G3   = zeros(h,w,d);
G31  = zeros(h,w,d);
G32  = zeros(h,w,d);
M1 = zeros(h,w,d);   % multiplier of (Yten - Xten - Sten);
M2 = zeros(h,w,d);   % multiplier of (nabla_3 Xten - G3);
M3 = zeros(h,w,d);   % multiplier of (nabla_1 G3 - G31);
M4 = zeros(h,w,d);   % multiplier of (nabla_2 G3 - G32);
M5 = zeros(h,w,d);   % multiplier of (Lten - Xten);
iter = 0;
tic
while iter<maxIter
    iter = iter +1;
    %% update S
    Sten = softthre(Yten-Xten+M1/mu, lambda/mu); % lambda need to fine tuned
    %% update X
    diffT_p  = diff_zT(mu*G3-M2,sizeD);
    numer1   = reshape( diffT_p + mu*(Yten(:)-Sten(:)) + M1(:) + mu*Lten(:) + M5(:), sizeD);
    x        = real( ifftn( fftn(numer1) ./ (mu*Eny_z + 2*mu ) ) );
    Xten     = reshape(x,[h,w,d]);
    %% update L
    Lten     =     proxF_tSVD_1(Xten-M5/mu,beta/mu,[]);
    %% update G3
    diffT_p  = diff_xT(mu*G31-M3,sizeD)+diff_yT(mu*G32-M4,sizeD);
    numer1   = reshape( diffT_p + mu*diff_z(Xten,sizeD)+M2(:), sizeD);
    x        = real( ifftn( fftn(numer1)./(mu*determ + mu ) ) );
    G3       = reshape(x,[h,w,d]);
    %% update G31 and G32
    G31      = softthre(reshape(diff_x(G3,sizeD),sizeD)+M3/mu,1/mu); % beta need to fine tuned
    G32      = softthre(reshape(diff_y(G3,sizeD),sizeD)+M4/mu,1/mu); 
    %% update multipler
    leq1 = Yten - Xten - Sten;
    leq2 = reshape(diff_z(Xten,sizeD),sizeD) - G3;
    leq3 = reshape(diff_x(G3,sizeD),sizeD) - G31;
    leq4 = reshape(diff_y(G3,sizeD),sizeD) - G32;
    leq5 = Lten - Xten;
    stopC1 = norm(leq1(:),'fro')/normD;
    stopC2 = norm(leq2(:),'fro')/normD;
    stopC4 = norm(leq4(:),'fro')/normD;
    if mod(iter,10)==0
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e')  ...
                ',Y-X-E=' num2str(stopC1,'%2.3e') ',||DX-X1||=' num2str(stopC2,'%2.3e')...
                ',|DZ-X3|' num2str(stopC4,'%2.3e')]);
    end
    if stopC1<tol && stopC2<tol
        break;
    else 
        M1 = M1 + mu*leq1;
        M2 = M2 + mu*leq2;
        M3 = M3 + mu*leq3;
        M4 = M4 + mu*leq4;
        M5 = M5 + mu*leq5;
        mu = min(mu*rho,1e6);
    end
end
