

%
% Moisture model simulation according to paper Kochanski et al., 2012.
% Code structure conforms as much as possible to 
%
%

%
%  This function will run the moisture model for one time step.
%  Synopsis: 
%            m = moisture_model(T, Q, P, m, w, r, dt)
%
%  Arguments:
%
%            T - the temperature in Kelvin
%            Q - the water vapor content (fraction, dimensionless)
%            P - the current surface atmospheric pressure (Pascals)
%            m - the current moisture content in the fuel (fraction,
%                          dimensionless)
%            r - the current rainfall intensity (mm/h, valid for current
%                          time step)
%            dt - the integration step [s]
%
%
%  Returns:
%
%            m - the new moisture values for the time t+dt
%

function m = moisture_model(T, Q, P, m, r, dt)

    r0 = 0.05;      % threshold rainfall [mm/h]
    rk = 8;         % saturation rain intensity [mm/h]
    Trk = 14 * 3600;% time constant for wetting model [s]
    S = 2.5;        % saturation intensity [dimensionless]
    
    Tk = [1, 10, 100]' * 3600;  % time lags for fuel classes [s]

    % saturated vapor pressure
    Pws = exp(54.842763 - 6763.22/T - 4.210 * log(T) + 0.000367*T + ...
          tanh(0.0415*(T - 218.8)) * (53.878 - 1331.22/T - 9.44523 * ...
          log(T) + 0.014025*T));
      
    % water vapor pressure
    Pw = P * Q / (0.622 + (1 - 0.622) * Q);
    
    % relative humidity (percent)
    H = 100 * Pw / Pws;
    
    % drying/wetting fuel equilibrium moisture contents
    Ed = 0.924*H^0.679 + 0.000499*exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - exp(-0.115*H));
    Ew = 0.618*H^0.753 + 0.000454*exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - exp(-0.115*H));
    
    Ed = Ed * 0.01;
    Ew = Ew * 0.01;
    
    % if rainfall is above threshold, apply saturation model
    if(r > r0)
    
        % equilibrium is equal to the saturation level
        equi = ones(size(m)) * S;
        
        % rlag is modified by the rainfall intensity
        rlag = ones(size(m)) * 1.0 / Trk * (1 - exp(- (r - r0) / rk));
    
    else
        
        % equilibrium is selected according to current moisture level
        equi = (m > Ed) * Ed + (m < Ew) * Ew + (m <= Ed) .* (m >= Ew) .* m;
        
        % the inverted time lag is constant according to fuel category
        rlag = 1 ./ Tk;
        
    end
    
    % select appropriate integration method (small change -> big errors in
    % the standard solution)
    change = dt * rlag;
    for i=1:length(m)
        if(change(i) < 0.01)
            m(i) = m(i) + (equi(i) - m(i)) * (1 - exp(-change(i)));
        else
            m(i) = m(i) + (equi(i) - m(i)) * change(i) * (1 - 0.5 * change(i));
        end
    end
        
