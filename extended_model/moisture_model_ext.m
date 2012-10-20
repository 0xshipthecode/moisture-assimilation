

%
% Modified and extended moisture model simulation (based on model in paper
% Kochanski et al., 2012).  Code structure is based on moisture_model.m
%
%
%
%  This function will run the moisture model for one time step.
%  Synopsis: 
%            m_ext = moisture_model(T, Tk, Q, P, m_ext, f_info, r, dt)
%  Arguments:
%
%            T - the temperature in Kelvin [differs between spatial
%                     locations, vector n x 1]
%            Tk - the nominal time lags for each fuel class (seconds)
%            Q - the water vapor content (fraction, dimensionless) [differs
%                      between spatial locations, vector n x 1]
%            P - the current surface atmospheric pressure (Pascals)
%                      [differs between spatial locations, vector n x 1]
%            m_ext - the extended state of the new moisture model, which
%                     includes the fuel moisture (dimensionless fraction),
%                     changes in rates (1/Tk,1/Trk - all in Hz) and target
%                     values (Ed, Ew, S - all dimensionless fractions).
%                     Total length is n*k+k+3, where k is the no. of fuel
%                     classes and n is the number of spatial locations.
%            f_info - the fuel type and location index for each fuel
%                   modeled in m_ext (1st col = f_type, 2nd col = f_loc)
%            r - the current rainfall intensity (mm/h, valid for current
%                          time step) [differs between spatial locations,
%                          vector n x 1]
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

function [m_ext, model_ids] = moisture_model_ext(T, Tk, Q, P, m_ext, f_info, r, dt)

    % number of fuels modelled (not necessarily at the same spatial
    % location)
    f_type = f_info(:,1);
    f_loc = f_info(:,2);

    k = length(Tk);                 % number of fuel classes
    n = max(f_loc);                 % number of fuel locations 
    r0 = 0.05;                      % threshold rainfall [mm/h]
    rk = 8;                         % saturation rain intensity [mm/h]
    Trk = 14 * 3600;                % time constant for wetting model [s]
    S = 2.5;                        % saturation intensity [dimensionless]
    
    nk = n * k;
    
    % first, we break the state vector into components
    m = m_ext(1:nk);
    dlt_Tk = m_ext(nk+1:nk+k);
    dlt_E = m_ext(nk+k+1);
    dlt_S = m_ext(nk+k+2);
    dlt_Trk = m_ext(nk+k+3);

    % saturated vapor pressure (at each location, size n x 1)
    Pws = exp(54.842763 - 6763.22./T - 4.210 * log(T) + 0.000367*T + ...
          tanh(0.0415*(T - 218.8)) .* (53.878 - 1331.22./T - 9.44523 * ...
          log(T) + 0.014025*T));
      
    % water vapor pressure (at each location, size n x 1)
    Pw = P .* Q ./ (0.622 + (1 - 0.622) * Q);
    
    % relative humidity (percent, at each location, size n x 1)
    H = 100 * Pw ./ Pws;
    
    % drying/wetting fuel equilibrium moisture contents (location specific,
    % n x 1)
    Ed = 0.924*H.^0.679 + 0.000499*exp(0.1*H) + 0.18*(21.1 + 273.15 - T).*(1 - exp(-0.115*H));
    Ew = 0.618*H.^0.753 + 0.000454*exp(0.1*H) + 0.18*(21.1 + 273.15 - T).*(1 - exp(-0.115*H));
    
    % rescale to range <0,1>, additionally add assimilated difference,
    % which is shared across spatial locations
    Ed = Ed * 0.01 + dlt_E;
    Ew = Ew * 0.01 + dlt_E;
    
    % where rainfall is above threshold (spatially different), apply
    % saturation model, equi and rlag are specific to fuel type and
    % location
    equi = m;  % copy over current equilibrium levels
    rlag = zeros(nk,1);
    model_ids = zeros(nk,1);
    
    % equilibrium is equal to the saturation level (assimilated)
    rain_model = r(f_loc) > r0;
    equi(rain_model) = (S + dlt_S);
    model_ids(rain_model) = 3;

    % rlag is modified by the rainfall intensity (assimilated)
    rlag(rain_model) = 1.0 / (Trk + dlt_Trk) .* (1 - exp(- (r(f_loc(rain_model)) - r0) / rk));
    
    % equilibrium is selected according to current moisture level
    equi(~rain_model & equi > Ed(f_loc)) = Ed(f_loc(~rain_model & equi > Ed(f_loc)));
    equi(~rain_model & equi < Ew(f_loc)) = Ew(f_loc(~rain_model & equi < Ew(f_loc)));
    model_ids(~rain_model) = (m(~rain_model) > Ed(f_loc(~rain_model))) * 1 ...
                           + (m(~rain_model) < Ew(f_loc(~rain_model))) * 2;
    model_ids(~model_ids) = 4;
        
    % the inverted time lag is constant according to fuel category
    % modified model: time constants for fuels are affected by
    % assimilation
    rlag(~rain_model) = 1 ./ (Tk(f_type(~rain_model)) + dlt_Tk(f_type(~rain_model)));        
    
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

