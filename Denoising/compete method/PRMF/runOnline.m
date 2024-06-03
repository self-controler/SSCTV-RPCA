load hall1-200;
matlabpool open;
X = normalize(XO);
[P Q L] = onlineRPMF(X, 2, 1, 1, 1e-2, ones(size(X)));
show(X, L, abs(X - L), [144 176]);
matlabpool close;