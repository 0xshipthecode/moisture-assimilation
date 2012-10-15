

%
%  This script runs the extended moisture model for one grid point with the
%  Extended Kalman filter.
%


% simulation time in hours (integration step for moisture is 1h)
t = (0:1:150)';
N = length(t);

% parameters of the simulation
T = 300;            % surface temperature, Kelvin
q = 0.005;          % water vapor content (dimensionless)
p = 101325;         % surface pressure, Pascals
n_k = 2;            % number of fuel categories
Tk = [10, 100]' * 3600;  % time lags for fuel categories
Ndim = 2*n_k + 3;

% external driving: rainfall characteristics
r = zeros(N,1);
r((t > 5) .* (t < 65) > 0) = 1.1; % 2mm of rainfall form hour 5 to 65

% measured moisture at given times
obs_time = [2, 20, 50, 110, 140]';   % observation time in hours
obs_moisture = [ 0.05,  0.045; ...
                 0.3, 0.2; ...
                 0.7, 0.6; ...
                 0.03, 0.04; ...
                 0.032, 0.035]; % measurements for the n_k fuel classes
N_obs = length(obs_time);
current_obs = 1;

% initialization of the Kalman filter
m_ext = zeros(Ndim,1);
m_ext(1:n_k) = 0.03;

P = eye(Ndim) * 0.01;   % error covariance of the initial guess

% Kalman filter Q (model error covariance) and R (measurement error covar)
Q = eye(Ndim) * 0.01;
R = eye(n_k) * 0.5;

% the observation operator is a n_k x Ndim matrix with I_(n_k) on the left
H = zeros(n_k, Ndim);
H(1:n_k,1:n_k) = eye(n_k);

% storage space for results (with filtering)
m_f = zeros(N, Ndim);
m_f(1, :) = m_ext';

m_n = zeros(N, Ndim); % (without filtering)
m_n(1, :) = m_ext';

% indicator of moisture model that is switched on at each time point
model_ids = zeros(N, n_k);

% storage for matrix traces
trP = zeros(N, 1);
trK = zeros(N, 1);
trS = zeros(N, 1);
trJ = zeros(N, 1);
sP = zeros(N, Ndim, Ndim);
sJ = zeros(N, Ndim, Ndim);

% predict & update loop
for i=2:N
    
    % compute the integration time step
    dt = (t(i) - t(i-1)) * 3600;
    
    % compute & store results for system without Kalman filtering
    m_n(i, :) = moisture_model_ext(T, Tk, q, p, m_n(i-1,:)', r(i), dt);

    % KALMAN PREDICT STEP
    
    % estimate new moisture mean based on last best guess (m)
    [m_pred, model_ids(i,:)] = moisture_model_ext(T, Tk, q, p, m_f(i-1,:)', r(i), dt);
    
    % update covariance matrix using the tangent linear model
    Jm = moisture_tangent_model_ext(T, Tk, q, p, m_f(i-1,:)', r(i), dt);
    sJ(i, :, :) = Jm;
    P = Jm*P*Jm' + Q;
    sP(i, :, :) = P;
    %trP(i) = trace(P);
    trP(i) = abs(prod(eig(P)));
    trJ(i) = prod(eig(Jm));
    
    % KALMAN UPDATE STEP (if measurement available) 
    if((current_obs <= N_obs) && (t(i) == obs_time(current_obs)))
        
        % acquire current measurement & move to next one
        m_measured = obs_moisture(current_obs, :)';
        current_obs = current_obs + 1;

        % innovation covariance (H=I due to direct observation)
        S = H*P*H' + R;
%        trS(i) = trace(S);
        trS(i) = prod(eig(S));
        
        % Kalman gain is inv(S) * P for this case (direct observation)
        K = P * H' / S;
        %trK(i) = trace(K);
        trK(i) = prod(eig(K(1:n_k,1:n_k)));
        
        % update step of Kalman filter to shift model state
        m_f(i,:) = m_pred + K*(m_measured - H*m_pred);
        
        % state error covariance is reduced by the observation
        P = P - K*S*K';
    
    else
        
        % if no observation is available, store the predicted value
        m_f(i,:) = m_pred;
        
    end
        
end

figure;
subplot(311);
plot(t, m_f(:,1), 'g-', 'linewidth', 2);
hold on;
plot(repmat(t, 1, 2), [m_f(:,1) - sqrt(sP(:, 1, 1)), m_f(:,1) + sqrt(sP(:, 1, 1))], 'gx');
plot(t, m_n(:,1), 'r-', 'linewidth', 2);
plot(t, r, 'k--', 'linewidth', 2);
plot(obs_time, obs_moisture(:,1), 'ko', 'markersize', 8, 'markerfacecolor', 'b');
h = legend('system + EKF', '$$m - \sigma$$', '$$m + \sigma$$', 'raw system', 'rainfall [mm/h]', 'observations');
set(h, 'interpreter', 'latex');
title('Plot of the evolution of the moisture model [EKF]', 'fontsize', 16);

% select time indices corresponding to observation times
[I,J] = ind2sub([N_obs, N], find(repmat(t', N_obs, 1) == repmat(obs_time, 1, N)));
subplot(312);
plot(t, log10(trP), 'b-', 'linewidth', 2);
hold on;
plot(t, log10(trJ), 'k-', 'linewidth', 2);
plot(obs_time, log10(trS(J)), 'ko', 'markerfacecolor', 'green', 'markersize', 6);
plot(obs_time, log10(trK(J)), 'ko', 'markerfacecolor', 'red', 'markersize', 6);
hold off;
legend('State', 'Jacobian', 'Innovation', 'Kalman gain');
title('Kalman filter: log(generalized variance) of covar/Kalman matrices vs. time [EKF]', 'fontsize', 16);

subplot(313);
plot(t, sP(:, 1, 1), 'r-', 'linewidth', 2);
hold on
plot(t, sP(:, 1, 2), 'g-', 'linewidth', 2);
plot(t, sP(:, 1, 3), 'b-', 'linewidth', 2);
plot(t, sP(:, 1, 4), 'k-', 'linewidth', 2);
plot(t, sP(:, 1, 5), 'm-', 'linewidth', 2);
hold off
legend('var(m)', 'cov(m,dT)', 'cov(m,dE)', 'cov(m,dS)', 'cov(m,dTr)');
title('Covariance between moisture and system parameters [EKF]');

figure;
subplot(311);
plot(repmat(t, 1, 2), m_f(:,[2,5]), 'linewidth', 2);
title('Time constant changes');
legend('dTk1', 'dTrk');
subplot(312);
plot(repmat(t, 1, 2), m_f(:, [3,4]), 'linewidth', 2);
legend('dE', 'dS');
title('Equilibrium changes');
subplot(313);
plot(t, model_ids, 'or');
title('Active submodel of the moisture model [EKF]');
