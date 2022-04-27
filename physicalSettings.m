function [c, rho, CST] = physicalSettings( T )
%  PHYSICALSETTINGS: Determine various physical settings for a given
%  temperature.
%
% [C, RHO, CST] = PHYSICALSETTINGS( T ) returns the speed of sound C, the
% density of air RHO, and a constant used for losses inside acoustic
% structures based on a temperature in degrees Celsius.
%
% References:
%
% 1. D. Keefe, "Acoustical wave propagation in cylindrical ducts:
%      Transmission line parameter approximations for isothermal and
%      nonisothermal boundary conditions," Journal of the Acoustical
%      Society of America, Vol. 75, No. 1, pp. 58-62, 1984.
%
% 2. A. Lefebvre, G. Scavone, J. Kergomard, "External Tonehole Interactions
%      in Woodwind Instruments," Acta Acustica united with Acustica,
%      Vol. 99, pp. 975-985, 2013.
%
% by Gary P. Scavone, McGill University, 2013-2021.

deltaT = T - 26.85;
c = 347.23 * ( 1 + 0.00166 * deltaT );        % speed of sound in air (m/s)
rho = 1.1769 * ( 1 - 0.00335 * deltaT );      % density of air (kg/m^3)
mu = 1.846*10^(-5) * ( 1 + 0.0025 * deltaT ); % shear viscosity coefficient (kg/m s)
gamma = 1.4017 * ( 1 - 0.00002 * deltaT );    % ratio of specific heats
sqrtPr = 0.8410 * ( 1 - 0.00002 * deltaT );   % Prandtl number (around 0.71)
CST = sqrt(mu/(rho*c)/2)*(1+((gamma-1)/sqrtPr)); % loss constant (see Ref. 2)