

%
% Modified and extended moisture model simulation (based on model in paper
% Kochanski et al., 2012).  Code structure is based on moisture_model.m
%
%
%
%  This function will run the moisture model for one time step.
%  Synopsis: 
%            m_ext = moisture_model(T, Q, P, m_ext, w, r, dt)
%
%  Arguments:
%
%            T - the temperature in Kelvin
%            Q - the water vapor content (fraction, dimensionless)
%            P - the current surface atmospheric pressure (Pascals)
%            m_ext - the extended state of the new moisture model, which
%                     includes the fuel moisture (dimensionless fraction),
%                     changes in rates (1/Tk,1/Trk - all in Hz) and target
%                     values (Ed, Ew, S - all dimensionless fractions).
%                     Total length is 2*k+3, where k is the no. of fuel
%                     classes.
%            r - the current rainfall intensity (mm/h, valid for current
%                          time step)
%            dt - the integration step [s]
%
%
%  Returns:
%
%            m_ext - the new moisture values (and state of extended model)
%                    for the time t+dt
%
%            model_id - an integer identifier of the model that is active
%                      (1 - drying, 2 - wetting, 3 - rain, 4 - dead zone)
%

function [m_ext, model_id] = moisture_model_ext(T, Q, P, m_ext, r, dt)

    k = (length(m_ext) - 4) / 2;    % number of fuel components
    r0 = 0.05;                      % threshold rainfall [mm/h]
    rk = 8;                         % saturation rain intensity [mm/h]
    Trk = 14 * 3600;                % time constant for wetting model [s]
    S = 2.5;                        % saturation intensity [dimensionless]
    
    Tk = [1, 10, 100]' * 3600;  % time lags for fuel classes [s]
    
    % first, we break the state vector into components
    m = m_ext(1:k);
    dlt_Tk = m_ext(k+1:2*k);
    dlt_Ed = m_ext(2*k+1);
    dlt_Ew = m_ext(2*k+2);
    dlt_S = m_ext(2*k+3);
    dlt_Trk = m_ext(2*k+4);

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
    
    % rescale to fractions
    % modification: extended model equilibria affected by assimilation
    Ed = Ed * 0.01 + dlt_Ed;
    Ew = Ew * 0.01 + dlt_Ew;
    
    % if rainfall is above threshold, apply saturation model
    if(r > r0)
    
        % equilibrium is equal to the saturation level
        % modification in extended model: saturation perturbed by
        % assimilation
        equi = ones(size(m)) * (S + dlt_S);
        model_id = 3 * ones(k, 1);
        
        % rlag is modified by the rainfall intensity
        % modification in extended model: wetting rate perturbed by
        % assimilation
        rlag = ones(size(m)) * 1.0 / (Trk + dlt_Trk) * (1 - exp(- (r - r0) / rk));
    
    else
        
        % equilibrium is selected according to current moisture level
        % note: equilibrium is modified by data assimilation above
        equi = m(1:k);
        equi(equi > Ed) = Ed;
        equi(equi < Ew) = Ew;
        model_id = (m(1:k) > Ed) * 1 + (m(1:k) < Ew) * 2 + (m(1:k) >= Ew) .* (m(1:k) <= Ed) * 4;
        
        % the inverted time lag is constant according to fuel category
        % modified model: time constants for fuels are affected by
        % assimilation
        rlag = 1 ./ (Tk + dlt_Tk);
        
    end
    
    % select appropriate integration method (small change -> big errors in
    % the standard solution)
    change = dt * rlag;
    for i=1:k
        if(change(i) < 0.01)
            m_ext(i) = m(i) + (equi(i) - m(i)) * (1 - exp(-change(i)));
        else
            m_ext(i) = m(i) + (equi(i) - m(i)) * change(i) * (1 - 0.5 * change(i));
        end
    end
    
    % modified model: delta variables do not evolve in model, only during
    % assimilation

