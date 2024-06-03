load hall1-200;
% matlabpool open;
X = normalize(XO);
[P,Q] = RPMF(X, 2, 1, 1, 1e-2);
show(X, P * Q, abs(X - P * Q), [144 176]);
% matlabpool close;