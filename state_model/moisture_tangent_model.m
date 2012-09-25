

%
% Tangent linear model to the moisture_model.m nonlinear moisture model.
%
%

%
%  This function will run the tangent linear moisture model for one time step.
%
%  Synopsis: 
%            m = moisture_model_tangent(T, Q, P, m, w, r, dt)
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
%            m - the Jacobian of the nonlinear model at previous (m)
%

function Jm = moisture_tangent_model(T, Q, P, m, r, dt)

    r0 = 0.05;      % threshold rainfall [mm/h]
    rk = 8;         % saturation rain intensity [mm/h]
    Trk = 14 * 3600;% time constant for wetting model [s]
    S = 2.5;        % saturation intensity [dimensionless]
    
    Tk = [1, 10, 100]' * 3600;  % time lags for fuel classes [s]
    
    % if rainfall is above threshold, apply saturation model
    if(r > r0)
    
        % rlag is modified by the rainfall intensity
        rlag = ones(size(m)) * 1.0 / Trk * (1 - exp(- (r - r0) / rk));
    
    else
        
        % the inverted time lag is constant according to fuel category
        rlag = 1 ./ Tk;
        
    end
    
    % select appropriate integration method (small change -> big errors in
    % the standard solution)
    change = dt * rlag;
    Jm = zeros(3,3);
    for i=1:3
        if(change(i) < 0.01)
            Jm(i,i) = exp(-change(i));
        else
            Jm(i,i) = 1.0 - change(i) * (1 - 0.5 * change(i));
        end
    end
        
