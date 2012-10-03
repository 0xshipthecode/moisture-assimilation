

%
% Tangent linear model to the moisture_model_ext.m extended nonlinear 
% moisture model.
%
%
%
%  This function will run the tangent linear moisture model for one time step.
%
%  Synopsis: 
%            m_ext = moisture_model_ext_tangent(T, Q, P, m_ext, w, r, dt)
%
%  Arguments:
%
%            T - the temperature in Kelvin
%            Q - the water vapor content (fraction, dimensionless)
%            P - the current surface atmospheric pressure (Pascals)
%            m - the extended state (see moisture_model_ext.m) including
%                current moisture content in the fuel (dimless. fraction)
%            r - the current rainfall intensity (mm/h, valid for current
%                          time step)
%            dt - the integration step [s]
%
%
%  Returns:
%
%            Jm_ext - the Jacobian of the nonlinear model at previous (m_ext)
%

function Jm_ext = moisture_tangent_model_ext(T, Q, P, m_ext, r, dt)

    k = (length(m_ext) - 3)/2;      % number of fuel components
    r0 = 0.05;                      % threshold rainfall [mm/h]
    rk = 8;                         % saturation rain intensity [mm/h]
    Trk = 14 * 3600;                % time constant for wetting model [s]
    S = 2.5;                        % saturation intensity [dimensionless]
    
    Tk = [1, 10, 100]' * 3600;  % time lags for fuel classes [s]
    
    % first, we break the state vector into components
    m = m_ext(1:k);
    dlt_Tk = m_ext(k+1:2*k);
    dlt_E = m_ext(2*k+1);
    dlt_S = m_ext(2*k+2);
    dlt_Trk = m_ext(2*k+3);
    
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
    Ed = Ed * 0.01 + dlt_E;
    Ew = Ew * 0.01 + dlt_E;
    
    % if rainfall is above threshold, apply saturation model
    if(r > r0)
    
        % equilibrium is equal to the saturation level
        equi = ones(size(m)) * (S + dlt_S);

        % rlag is modified by the rainfall intensity
        rlag = ones(size(m)) * 1.0 / (Trk + dlt_Trk) * (1 - exp(- (r - r0) / rk));
    
    else
        
        % the equilibrium level is given depending on the current moisture
        % state and is Ed (drying equilibrium) if moisture is above Ed, Ew
        % (wetting equilibrium) if model is below Ew or the moisture value
        % itself if in between).
        equi = m;
        equi(m > Ed) = Ed;
        equi(m < Ew) = Ew;

        % the inverted time lag is constant according to fuel category
        rlag = 1 ./ (Tk + dlt_Tk);

    end
    
    % select appropriate integration method (small change -> big errors in
    % the standard solution)
    change = dt * rlag;
    Jm_ext = zeros(2*k+3);
    for i=1:k
        
        if(change(i) < 0.01)
            
            % partial m_i/partial m_i
            Jm_ext(i,i) = exp(-change(i));
            
            % precompute partial m_i/partial change
            dmi_dchng = (equi(i) - m(i)) * exp(-change(i));
            
            % precompute partial m_i/partial equi
            dmi_dequi = (1.0 - exp(-change(i)));

        else
            
            % partial dm_i/partial m_i
            Jm_ext(i,i) = 1.0 - change(i) * (1 - 0.5 * change(i));
            
            % partial m_i/partial change
            dmi_dchng = (equi(i) - m(i)) * (1 - change(i));
            
            % partial m_i/partial equi
            dmi_dequi = change(i) * (1 - 0.5 * change(i));
                        
        end
        
        
        % branch according to the currently active model
        if(r <= r0)
            
            % drying/wetting model active

            % partial m_i/partial delta_Tk
            Jm_ext(i,k+i) = dmi_dchng * (-dt) * (Tk(k) + dlt_Tk(k))^(-2);

            % if drying/wetting model active, jacobian entry is nonzero;
            % it is zero if the 'dead zone' model is active
            if((m(i) > Ed) || (m(i) < Ew))
                Jm_ext(i,2*k+1) = dmi_dequi;
            end

        else

            % rain model active

            % partial m_i/partial deltaS
            Jm_ext(i,2*k+2) = dmi_dequi;

            % partial m_i/partial deltaTkr
            Jm_ext(i,2*k+3) = dmi_dchng * dt * (exp(-(r - r0)/rk) - 1) * (Trk + dlt_Trk)^(-2);

        end
        
        % delta_Tk for each fuel have no dependencies except previous delta_Tk
        Jm_ext(k+i, k+i) = 1.0;        
        
    end
    
    % the equilibrium constants 
    Jm_ext(2*k+1,2*k+1) = 1.0;
    Jm_ext(2*k+2,2*k+2) = 1.0;
    Jm_ext(2*k+3,2*k+3) = 1.0;
