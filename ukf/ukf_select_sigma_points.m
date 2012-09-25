
% This function selects a sigma ensemble that conforms to the requirements
% for deterministic sampling in the UKF filter.  This method automatically
% generates 2*N points (in addition to the state vector itself), where N is
% the dimension of the state vector.
%
%  synopsis: [X,W] = ukf_select_sigma_points(x, Sigma, W0)
%
%   x - the current state vector
%   Sigma - the current covariance matrix of the state
%   W0 - the weight of the state vector in the sigma point set
%

function [X,W] = ukf_select_sigma_points(x, Sigma, W0)

    N = size(x,1);
    
    X = zeros(N, N);
    W = zeros(2*N+1,1);
    
    W(:) = (1 - W0)/(2*N);
    W(2*N+1) = W0;
    X(:, 2*N+1) = x;
    
    % matlab decomposition is A'*A, so we use the rows
    SF = chol(Sigma);
    for i=1:N
       X(:,i) = x + sqrt(0.5/W(i)) * SF(i,:)';
       X(:,i+N) = x - sqrt(0.5/W(i)) * SF(i,:)';
    end
