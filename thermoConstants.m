function [c, rho, gamma, lv, Pr] = thermoConstants( T )
% THERMOCONSTANTS: Determine various thermodynamic constants for a
% specified temperature.
%
% [C, RHO, GAMMA, LV, PR] = THERMOCONSTANTS( T ) returns the speed of sound
% C, the density of air RHO, the ratio of specific heats GAMMA, the viscous
% characteristic length LV, and the Prandtl number PR, based on a
% temperature T in degrees Celsius (default T = 20 degrees Celsius).
%
% Reference:
%
% 1. D. Keefe, "Acoustical wave propagation in cylindrical ducts:
%      Transmission line parameter approximations for isothermal and
%      nonisothermal boundary conditions," Journal of the Acoustical
%      Society of America, Vol. 75, No. 1, pp. 58-62, 1984.
%
% by Gary P. Scavone, McGill University, 2013-2024.

if ~exist( 'T', 'var')
  T = 20;
end

deltaT = T - 26.85;
c = 347.23 * ( 1 + 0.00166 * deltaT );        % speed of sound in air (m/s)
rho = 1.1769 * ( 1 - 0.00335 * deltaT );      % density of air (kg/m^3)
mu = 1.846*10^(-5) * ( 1 + 0.0025 * deltaT ); % shear viscosity coefficient (kg/m s)
gamma = 1.4017 * ( 1 - 0.00002 * deltaT );    % ratio of specific heats
Pr = (0.8410 * ( 1 - 0.00002 * deltaT ))^2;   % sqrt of Prandtl number (around 0.71)
lv = mu/(rho*c);
