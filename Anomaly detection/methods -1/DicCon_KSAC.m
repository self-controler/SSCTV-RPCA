function Dictionary = DicCon_KSAC(X, K, P)
%
%   Dictionary construction for LRASR 
%
% Inputs:
%   X  - 2D data matrix (p x N) 
%   K - number of clusters
%   P - number of selected pixels
% Outputs:
%   Dictionary - Detector output (1 x N)
%
%Author: Yang Xu

% K-means 
[K1 K2]=kmeans(X,K);
Dictionary = [];
for i=1:K

    st1=find(K1==i);
    
    if length(st1)<P
        continue;
    end
    temp=X(:,st1);
    %kr=RxDetector(X(:,st1));
    kr=KSAC(X(:,st1));
    %[d1 d2]=sort(kr,'ascend');
    [d1 d2]=sort(kr,'descend');
    Dictionary=[Dictionary ,temp(:,d2(1:P))];
end
    