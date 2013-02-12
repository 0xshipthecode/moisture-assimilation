
%  This function computes the AAR clustering from the connectivity matrix
%  C.  The number of clusters returned is controlled by K.  The clustering
%  method returns disjoint clusters but the clusters need not cover the
%  entire set of voxels.
%
%  synopsis: C = cluster_arr(CC, K)
%
%   CC - input connectivity matrix (must be symmetric) (N x N)
%   K - the number of clusters to return
%   C - a matrix of N x K elements
%
%

function C = cluster_aar(CC, K)

    % extract eigenvalues and eigenvectors
    [V, D] = eig(CC);
    D = diag(D);
    
    % sort eigenvalues in descending order
    [~, ord] = sort(D, 'descend');
    V = V(:, ord);
    
    % extract first K eigvecs
    X = V(:, 1:K);
    nzs = sum(X.^2, 2) > 0;
    
    [~, R] = rotatefactors(X(nzs, :), 'Method', 'varimax', ...
             'Maxit', 100000, 'Normalize', 'off');

    % rotate the eigenvectors to the varimax-determined base
    Y = X * R;
    
    % find the inliers under this hypothesis
    C = nonmax_suppress_builtin(abs(Y));
    for i=1:K
        
        Yi = abs(Y(:,i)) .* C(:,i);
        Yi = Yi ./ sqrt(sum(Yi.^2));

        % fit the discrete indicator vector 
        C(:,i) = fit_best_indicator_fast_builtin(Yi);
        
    end

    
%
% This function performs nonmax_suppression on the matrix X (row-wise),
% where the maxima are taken from the absolute values of X.
%
% synopsis: Xs = nonmax_suppress(X)
%
%
function Xs = nonmax_suppress_builtin(X)
    
    [N,K] = size(X);
    
    [~, Ypos] = sort(X, 2);
    Xs = zeros(N,K);
    
    for i=1:N
        Xs(i,Ypos(i,K)) = 1;
    end
        
    
function z = fit_best_indicator_fast_builtin(w)

    w_backup = w;
    
    w = sort(w);
    w = w(end:-1:1);
    
    maxL = find(w > 0, 1, 'last');
    
    if(isempty(maxL))
        z = zeros(size(w));
        return
    end
    
    dists = ((1:maxL) + 1) ./ (1:maxL) - 2 ./ sqrt(1:maxL) .* cumsum(w(1:maxL))';
    [~, ThrNdx] = min(dists);
    thr = w(ThrNdx);
    
    % select the elements which belong to the cluster
    z = zeros(size(w));
    z(w_backup >= thr) = 1;
    