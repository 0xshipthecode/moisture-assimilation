

%
%  This script runs the moisture model for one grid point.
%



% time in hours (integration step for moisture is 1h)
t = (0:0.1:400)';
N = length(t);

% parameters of the simulation
T = 305; % Kelvin
Q = 0.005; % water vapor content (dimensionless)
P = 101325;  % Pascals
m = 0.04 * ones(3,1);
n_k = length(m);  % number of fuel categories


% prepare rainfall characteristics
r = zeros(N,1);
r((t > 5) .* (t < 105) > 0) = 1.0; % 1mm of rainfall during 

m_t = zeros(N, n_k);
m_t(1, :) = m';
for i=2:N
    dt = (t(i) - t(i-1)) * 3600;
    m_new = moisture_model(T, Q, P, m, r(i), dt);
    m_t(i, :) = m_new;
    m = m_new;
end

plot(repmat(t, 1, n_k), m_t);