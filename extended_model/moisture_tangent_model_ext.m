

%
% Tangent linear model to the moisture_model_ext.m extended nonlinear 
% moisture model.
%
%
%
%  This function will run the tangent linear moisture model for one time step.
%
%  Synopsis: 
%            m_ext = moisture_model_ext_tangent(T, Tk, Q, P, m_ext, f_info, w, r, dt)
%
%  Arguments:
%
%            T - the temperature in Kelvin
%            Tk - the nominal time lags for each fuel class (seconds)
%            Q - the water vapor content (fraction, dimensionless)
%            P - the current surface atmospheric pressure (Pascals)
%            m - the extended state (see moisture_model_ext.m) including
%                current moisture content in the fuel (dimless. fraction)
%            f_info - the fuel type and location index for each fuel
%                   modeled in m_ext (1st col = f_type, 2nd col = f_loc)
%            r - the current rainfall intensity (mm/h, valid for current
%                          time step)
%            dt - the integration step [s]
%
%
%  Returns:
%
%            Jm_ext - the Jacobian of the nonlinear model at previous (m_ext)
%

function Jm_ext = moisture_tangent_model_ext(T, Tk, Q, P, m_ext, f_info, r, dt)

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
    
    % saturated vapor pressure
    Pws = exp(54.842763 - 6763.22./T - 4.210 * log(T) + 0.000367*T + ...
          tanh(0.0415*(T - 218.8)) .* (53.878 - 1331.22./T - 9.44523 * ...
          log(T) + 0.014025*T));
      
    % water vapor pressure
    Pw = P .* Q ./ (0.622 + (1 - 0.622) * Q);
    
    % relative humidity (percent)
    H = 100 * Pw ./ Pws;
    
    % drying/wetting fuel equilibrium moisture contents
    Ed = 0.924*H.^0.679 + 0.000499*exp(0.1*H) + 0.18*(21.1 + 273.15 - T).*(1 - exp(-0.115*H));
    Ew = 0.618*H.^0.753 + 0.000454*exp(0.1*H) + 0.18*(21.1 + 273.15 - T).*(1 - exp(-0.115*H));
    
    % rescale to fractions
    % modification: extended model equilibria affected by assimilation
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
    equi(rain_model) = S + dlt_S;

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
    Jm_ext = zeros(nk + k + 3);
    for i=1:nk
        
        if(change(i) < 0.01)
            
            % partial m_i/partial m_i
            if(model_ids(i) < 4)
                Jm_ext(i,i) = exp(-change(i));
            else
                Jm_ext(i,i) = 1.0;
            end
            
            % precompute partial m_i/partial change
            dmi_dchng = (equi(i) - m(i)) * exp(-change(i));
            
            % precompute partial m_i/partial equi
            dmi_dequi = (1.0 - exp(-change(i)));

        else
            
            % partial dm_i/partial m_i
            if(model_ids(i) < 4)
                Jm_ext(i,i) = 1.0 - change(i) * (1 - 0.5 * change(i));
            else
                Jm_ext(i,i) = 1.0;
            end
            
            % partial m_i/partial change
            dmi_dchng = (equi(i) - m(i)) * (1.0 - change(i));
            
            % partial m_i/partial equi
            dmi_dequi = change(i) * (1 - 0.5 * change(i));
                        
        end
        
        
        % branch according to the currently active model
        if(r <= r0)
            
            % drying/wetting model active
            if((m(i) > Ed(f_loc(i))) || (m(i) < Ew(f_loc(i))))

                % partial m_i/partial delta_Tk
                Jm_ext(i,nk+i) = dmi_dchng * (-dt) * (Tk(i) + dlt_Tk(i))^(-2);
                
                % partial m_i/partial delta_E
                Jm_ext(i,nk+k+1) = dmi_dequi;
            end

        else

            % rain model active

            % partial m_i/partial deltaS
            Jm_ext(i,nk+k+2) = dmi_dequi;

            % partial m_i/partial deltaTkr
            Jm_ext(i,nk+k+3) = dmi_dchng * dt * (exp(-(r - r0)/rk) - 1) * (Trk + dlt_Trk)^(-2);

        end
                
    end
    
    % delta_Tk for each fuel have no dependencies except previous delta_Tk
    for i=1:k
        Jm_ext(nk+i, nk+i) = 1.0;
%        Jm_ext(nk+i, nk+i) = 0.99;
    end
    
    % the equilibrium constants 
    Jm_ext(nk+k+1, nk+k+1) = 1.0;
    Jm_ext(nk+k+2, nk+k+2) = 1.0;
    Jm_ext(nk+k+3, nk+k+3) = 1.0;
%     Jm_ext(nk+k+1, nk+k+1) = 0.99;
%     Jm_ext(nk+k+2, nk+k+2) = 0.99;
%     Jm_ext(nk+k+3, nk+k+3) = 0.99;
