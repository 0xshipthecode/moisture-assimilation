



%
%  This script runs the moisture model for one grid point including the
%  Extended Kalman filter.
%
%


% simulation time in hours (integration step for moisture is 1h)
t = (0:0.1:200)';
N = length(t);

% parameters of the simulation
T = 300;            % surface temperature, Kelvin
q = 0.005;          % water vapor content (dimensionless)
p = 101325;         % surface pressure, Pascals
n_k = 3;            % number of fuel categories

% external driving: rainfall characteristics
r = zeros(N,1);
r((t > 5) .* (t < 65) > 0) = 1.1; % 2mm of rainfall form hour 5 to 65

% measured moisture at given times
obs_time = [2, 20, 50, 100]';   % observation time in hours
obs_moisture = [ 0.05, 0.05, 0.05; ...
                 0.3, 1.2, 1.1; ...
                 0.7, 1.4, 1.2; ...
                 0.03, 0.5, 0.7]; % measurements for the 3 fuel classes
N_obs = length(obs_time);
current_obs = 1;

% initialization of the Kalman filter
m = 0.03 * ones(3,1);       % initial mean value
P = eye(n_k, n_k) * 0.1;   % error covariance of the initial guess

% Kalman filter Q (model error covariance) and R (measurement error covar)
Q = eye(n_k, n_k) * 0.1;
R = eye(n_k, n_k) * 10.0;

% storage space for results (with filtering)
m_f = zeros(N, n_k);
m_f(1, :) = m';

m_n = zeros(N, n_k); % (without filtering)
m_n(1, :) = m;

% storage for matrix traces
trP = zeros(N, 1);
trK = zeros(N, 1);
trS = zeros(N, 1);

% predict & update loop
for i=2:N
    
    % compute the integration time step
    dt = (t(i) - t(i-1)) * 3600;
    
    % compute & store results for system without Kalman filtering
    m_n(i, :) = moisture_model(T, q, p, m_n(i-1,:)', r(i), dt);

    % KALMAN PREDICT STEP
    
    % estimate new moisture mean based on last best guess (m)
    m_pred = moisture_model(T, q, p, m_f(i-1,:)', r(i), dt);
    
    % update covariance matrix using the tangent linear model
    Jm = moisture_tangent_model(T, q, p, m_f(i-1,:)', r(i), dt);
    P = Jm*P*Jm' + Q;
    %trP(i) = trace(P);
    trP(i) = prod(eig(P));
    
    % KALMAN UPDATE STEP (if measurement available) 
    if((current_obs <= N_obs) && (t(i) == obs_time(current_obs)))
        
        % acquire current measurement & move to next one
        m_measured = obs_moisture(current_obs, :)';
        current_obs = current_obs + 1;

        % innovation covariance (H=I due to direct observation)
        S = P + R;
%        trS(i) = trace(S);
        trS(i) = prod(eig(S));
        
        % Kalman gain is inv(S) * P for this case (direct observation)
        K = inv(S) * P;
        %trK(i) = trace(K);
        trK(i) = prod(eig(K));
        
        % update step of Kalman filter to shift model state
        m_f(i,:) = m_pred + K*(m_measured - m_pred);
        
        % state error covariance is reduced by the observation
        P = P - K*S*K';
    
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
plot(obs_time, obs_moisture(:,1), 'ko', 'markersize', 8, 'markerfacecolor', 'b');
legend('system + EKF', 'raw system', 'rainfall [mm/h]', 'observations');
title('Plot of the evolution of the moisture model', 'fontsize', 16);

% select time indices corresponding to observation times
[I,J] = ind2sub([N_obs, N], find(repmat(t', N_obs, 1) == repmat(obs_time, 1, N)));
subplot(212);
plot(t, log10(trP), 'linewidth', 2);
hold on;
plot(obs_time, log10(trS(J)), 'ko', 'markerfacecolor', 'green', 'markersize', 10);
plot(obs_time, log10(trK(J)), 'ko', 'markerfacecolor', 'red', 'markersize', 10);
hold off;
legend('State', 'Innovation', 'Kalman gain');
title('Kalman filter: log(generalized variance) of covar/Kalman matrices vs. time', 'fontsize', 16);


