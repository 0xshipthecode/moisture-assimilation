
%
%
%  A check of the moisture model mapping.
%
%

T = 300;            % surface temperature, Kelvin
q = 0.01;           % water vapor content (dimensionless)
p = 101325;         % surface pressure, Pascals
n_k = 3;            % number of fuel categories

dt = 3600;

m_orig = [repmat( (1:1000)' * 0.001, 1, 3), zeros(1000, 7)];
m_map = zeros(1000,10);
model_id = zeros(1000, 3);
for i=1:1000
    [m_map(i, :), model_id(i,:)] = moisture_model_ext(T, q, p, m_orig(i,:)', 0, dt);
end

subplot(211);
plot(m_orig(:,1), m_map(:,1), '-');
subplot(212);
plot(m_orig(:,1), model_id(:,1), 'k-');