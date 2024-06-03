function [result] = KSAC2(M)
%RX anomaly detector
%   RxDetector performs the RX anomaly detector
%
% Usage
%   [result] = hyperRxDetector(M)
% Inputs
%   M  - 2D data matrix (p x N)
% Outputs
%   result - Detector output (1 x N)   ---->> 186*10000
%   sigma - Covariance matrix (p x p)
%   sigmaInv - Inverse of covariance matrix (p x p)

% Remove the data mean
[p, N] = size(M);
result = zeros(N, 1);
c=1e6;
for i=1:N
    %result(i) = sum(M(:,i).*Y)/sqrt(sum(M(:,i).^2)*sum(Y.^2));
    sumL=0;
    Y=M(:,i);
    for j=1:N
        if j~=i
            sumL=sumL+sum((Y-M(:,i)).^2);
        end
    end
    result(i) = exp(-sumL)/2/c;
end
result = abs(result);

return;