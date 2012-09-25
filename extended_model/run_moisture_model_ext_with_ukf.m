

%
%  This script runs the extended moisture model for one grid point with the
%  Unscented Kalman Filter.
%


% simulation time in hours (integration step for moisture is 1h)
t = (0:1:500)';
N = length(t);

% parameters of the simulation
T = 300;            % surface temperature, Kelvin
q = 0.005;          % water vapor content (dimensionless)
p = 101325;         % surface pressure, Pascals
n_k = 3;            % number of fuel categories

% external driving: rainfall characteristics
r = zeros(N,1);
r((t > 5) .* (t < 65) > 0) = 2; % 2mm of rainfall form hour 5 to 65

% measured moisture at given times
obs_time = [2, 20, 50, 100, 140]';   % observation time in hours
obs_moisture = [ 0.05, 0.05, 0.05; ...
                 0.3, 1.2, 1.1; ...
                 0.7, 1.4, 1.2; ...
                 0.03, 0.5, 0.7; ...
                 0.032, 0.4, 0.6]; % measurements for the 3 fuel classes
N_obs = length(obs_time);
current_obs = 1;

% initialization of the Kalman filter
m_ext = zeros(10,1);
m_ext(1:3) = 0.03 * ones(3,1);

P = eye(10) * 0.001;   % error covariance of the initial guess

% Kalman filter Q (model error covariance) and R (measurement error covar)
Q = eye(10) * 0;
R = eye(3) * 0.1;

% the observation operator is a 3 x 10 matrix with [I_3, 0]
H = zeros(3, 10);
H(1:3,1:3) = eye(3);

% storage space for results (with filtering)
m_f = zeros(N, 10);
m_f(1, :) = m_ext';

m_n = zeros(N, 10); % (without filtering)
m_n(1, :) = m_ext;

% indicator of moisture model that is switched on at each time point
model_ids = zeros(N, 3);

% storage for matrix traces
trP = zeros(N, 1);
trK = zeros(N, 1);
trS = zeros(N, 1);

% W0 is a UKF parameter affecting the sigma point distribution
W0 = -0.9;
Ndim = size(m_ext, 1);
Npts = Ndim * 2 + 1;

% Storage of the sigma points
sigma_pts = zeros(N, Ndim, Npts);

% predict & update loop
for i = 2:N
    
    % compute the integration time step
    dt = (t(i) - t(i-1)) * 3600;
    
    % draw the sigma points
    [m_sigma, w] = ukf_select_sigma_points(m_f(i-1,:)', P, W0);
    sigma_pts(i, :, :) = m_sigma;
    
    % compute & store results for system without Kalman filtering
    m_n(i, :) = moisture_model_ext(T, q, p, m_n(i-1,:)', r(i), dt);

    % UKF prediction step - run the sigma set through the nonlinear
    % function
    
    % estimate new moisture mean based on last best guess (m)
    m_sigma_1 = zeros(Ndim, Npts);
    for n=1:Npts
        m_sigma_1(:,n) = moisture_model_ext(T, q, p, m_sigma(:,n), r(i), dt);
    end
    [~, model_ids(i,:)] = moisture_model_ext(T, q, p, m_sigma(:,end), r(i), dt);
    
    % compute the prediction mean x_mean(i|i-1)
    m_pred = sum(m_sigma_1 * diag(w), 2);
    
    % estimate covariance matrix using the sigma point set
    sqrtP = (m_sigma_1 - repmat(m_pred, 1, Npts)) * diag(sqrt(w));
    P = Q + sqrtP * sqrtP';
    
    trP(i) = abs(prod(eig(P)));
    
    % KALMAN UPDATE STEP (if measurement available) 
    if((current_obs <= N_obs) && (t(i) == obs_time(current_obs)))
        
        % acquire current measurement & move to next one
        y_measured = obs_moisture(current_obs, :)';
        current_obs = current_obs + 1;

        % run the observation model on the sigma point ensemble
        % since the model is linear, I could actuall propagate m_pred
        Y = H * m_sigma_1;
        y_pred = sum(Y * diag(w), 2);
        
        % innovation covariance (H=I due to direct observation)
        sqrtS = (Y - repmat(y_pred, 1, Npts)) * diag(sqrt(w));
        S = sqrtS * sqrtS' + R; 
        trS(i) = prod(eig(S));
        
        % the cross covariance of state & observation errors
        Cxy = (m_sigma_1 - repmat(m_pred, 1, Npts)) * diag(w) * (Y - repmat(y_pred, 1, Npts))';
        
        % Kalman gain is inv(S) * P for this case (direct observation)
        K = Cxy * inv(S);
        %trK(i) = trace(K);
        trK(i) = prod(eig(K(1:3,1:3)));
        
        % update step of Kalman filter to shift model state
        m_f(i,:) = m_pred + K*(y_measured - y_pred);
        
        % state error covariance is reduced by the observation
        P = P - K*S*K';
    
        % replace the stored covariance by the updated covariance after
        % processing the measurement
        trP(i) = abs(prod(eig(P)));
        
    else
        
        % if no observation is available, store the predicted value
        m_f(i,:) = m_pred;
        
    end
        
end

figure;
subplot(211);
plot(t, m_f(:,1), 'g-', 'linewidth', 2);
hold on;
plot(t, m_n(:,1), 'r-', 'linewidth', 2);
plot(t, r, 'k--', 'linewidth', 2);
plot(t, model_ids(:,1), 'kx');
plot(obs_time, obs_moisture(:,1), 'ko', 'markersize', 8, 'markerfacecolor', 'b');
legend('system + UKF', 'raw system', 'rainfall [mm/h]', 'model id', 'observations');
title('Plot of the evolution of the moisture model', 'fontsize', 16);

% select time indices corresponding to observation times
[I,J] = ind2sub([N_obs, N], find(repmat(t', N_obs, 1) == repmat(obs_time, 1, N)));
subplot(212);
plot(t, log10(trP), 'linewidth', 1.5);
hold on;
%plot(t, log10(trJ), 'k-', 'linewidth', 1.5);
plot(obs_time, log10(trS(J)), 'ko', 'markerfacecolor', 'green', 'markersize', 10);
plot(obs_time, log10(trK(J)), 'ko', 'markerfacecolor', 'red', 'markersize', 10);
hold off;
legend('State', 'Innovation', 'Kalman gain');
title('Kalman filter: log(generalized variance) of covar/Kalman matrices vs. time', 'fontsize', 16);

figure;
subplot(311);
plot(repmat(t, 1, 4), m_f(:,[4:6,10]), 'linewidth', 1.5);
title('Time constant changes');
legend('dTk1', 'dTk2', 'dTk3', 'dS');
subplot(312);
plot(repmat(t, 1, 3), m_f(:, [7, 8, 9]), 'linewidth', 1.5);
legend('dEd', 'dEw', 'dS');
title('Equilibrium changes');
subplot(313);
plot(t, model_ids, 'or');
title('Active submodel of the moisture model');

figure;
rng = 40:200;
plot(t(rng), m_n(rng,1), 'k-', 'linewidth', 2);
hold on
plot(t(rng), m_f(rng,1), 'r-', 'linewidth', 2);
plot(repmat(t(rng), 1, 3), squeeze(sigma_pts(rng, 1, [1,2,11])), 'x');
hold off
title('Plot of the 1st (fast) fuel with sigma point projections');
legend('Model', 'Model+UKF', '\sigma-points');
