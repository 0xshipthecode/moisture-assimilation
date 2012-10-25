
%
%  Compute the equilibrium moisture field for many spatial locations
%  simultaneously given the ambient environmental conditions.
%
%  synopsis: [Ed,Ew] = equilibrium_moisture(P, Q, T)
%
%                      P - surface pressure [Pa]
%                      Q - water vapor ratio [-]
%                      T - surface temperature [K]
%                      Ed, Ew - drying/wetting equilibrium [dimensionless,
%                                                fraction between 0 and 1]
%
%  Note: shape of Ed, Ew will be the same as the shape of P,Q,T.  The
%  shapes of inputs P,Q,T must be equal.
%
%

function [Ed,Ew] = equilibrium_moisture(P, Q, T)

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

    % remap values
    Ed = Ed * 0.01;
    Ew = Ew * 0.01;