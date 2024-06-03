function [result, sigma, sigmaInv] = RxDetector(M)
%RX anomaly detector
%   RxDetector performs the RX anomaly detector
%
% Usage
%   [result] = hyperRxDetector(M)
% Inputs
%   M  - 2D data matrix (p x N)
% Outputs
%   result - Detector output (1 x N)
%   sigma - Covariance matrix (p x p)
%   sigmaInv - Inverse of covariance matrix (p x p)

% Remove the data mean
[p, N] = size(M);
mMean = mean(M, 2);
M = M - repmat(mMean, 1, N);

% Compute covariance matrix
sigma = hyperCov(M);
delta=1e-5;%regularizaton parameter
sigmaInv = inv(sigma+delta*eye(size(sigma)));

result = zeros(N, 1);
for i=1:N
    result(i) = M(:,i).'*sigmaInv*M(:,i);
end
result = abs(result);

return;