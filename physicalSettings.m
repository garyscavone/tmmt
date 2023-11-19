function [c, rho, gamma, lv, Pr, alphacm] = physicalSettings( T, f, humidity, P0 )
%  PHYSICALSETTINGS: Determine various physical variables and loss settings
%  for a specified temperature, frequencies, humidity and atmospheric
%  pressure.
%
% [C, RHO, WALLCST, ALPHACM] = PHYSICALSETTINGS( T, F, HUMIDITY, P0 )
% returns the speed of sound C, the density of air RHO, and a constant
% (WALLCST) used to account for wall losses inside acoustic structures
% based on a temperature in degrees Celsius. If the optional input argument
% F is included, a vector ALPHACM (of size equal to F) will be returned
% that can be used to include attenuation due to classical and molecular
% effects. If the temperature T is not specified, a default value of 20 C
% is assumed. If either the relative HUMIDITY (as percentage) or
% atmospheric pressure (P0, in Pascals) is not specified (needed to compute
% ALPHACM), values of 40% and 101325 Pascals will be used by default.
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
% 3. M. van Walstijn, M. Campbell, J. Kemp, D. Sharp, "Wideband
%      Measurements of the Acoustic Impedance of Tubular Objects," Acta
%      Acustica united with Acustica, Vol. 91, pp. 590-604, 2005.
%
% by Gary P. Scavone, McGill University, 2013-2022.

if ~exist( 'T', 'var')
  T = 20;
end

if nargin > 1 && ~isvector(f)
  error( 'f should be a 1D vector of frequency values in Hz.' );
end

deltaT = T - 26.85;
c = 347.23 * ( 1 + 0.00166 * deltaT );        % speed of sound in air (m/s)
rho = 1.1769 * ( 1 - 0.00335 * deltaT );      % density of air (kg/m^3)
mu = 1.846*10^(-5) * ( 1 + 0.0025 * deltaT ); % shear viscosity coefficient (kg/m s)
gamma = 1.4017 * ( 1 - 0.00002 * deltaT );    % ratio of specific heats
sqrtPr = 0.8410 * ( 1 - 0.00002 * deltaT );   % sqrt of Prandtl number (around 0.71)
%wallCst = sqrt(mu/(rho*c)/2)*(1+((gamma-1)/sqrtPr)); % wall loss constant (see Ref. 2)
lv = mu/(rho*c);
Pr = sqrtPr^2;
if nargin < 2
  alphacm = [];
  return;
end

if ~exist( 'humidity', 'var')
  humidity = 40;  % room humidity as percentage
end

if ~exist( 'P0', 'var')
  P0 = 101325;    % atmospheric pressure at sea level in Pascals
end

omega = 2 * pi * f;
muB = 0.6 * mu;  % bulk viscosity
alphac = omega.^2 * mu * (4/3 + muB/mu + (gamma-1)/sqrtPr^2) / (2*rho*c^3); % classical attenuation

pv = 0.065773*T^3 + 0.1445*T^2 + 59.34*T + 560.54; % vapor pressure of water
cv1oR = (-0.000002778*T^3 + 0.0007857*T^2 + 0.08599*T + 3.883)*10^-3; % ratio of specific-heat to gas constant for O2
cv2oR = (-0.00000009259*T^3 + 0.00035596*T^2 + 0.02212*T + 0.5525)*10^-3; % ratio of specific-heat to gas constant for N2
alphavgm = pi*(gamma-1)^2/(2*gamma);
alphav1g = alphavgm * cv1oR;  % maximum absorption per wavelength for O2
alphav2g = alphavgm * cv2oR;  % maximum absorption per wavelength for N2

h = 0.01 * humidity * pv / P0;  % humidity factor
G = 4.41*10^6*h*(0.05 + 100*h)/(0.391 + 100*h);
T1 = 293.16;
TA = 300 + deltaT; % absolute temperature in Kelvin
F = 6.142*((T1/TA)^(1/3) - 1);
tau1 = 1/(2*pi*(24 + G)); % relaxation time for oxygen in air (21%)
tau2 = sqrt(TA/T1)/(2*pi*(9 + (3.5*10^4)*h*exp(-F))); % relaxation time for nitrogen in air (78%)

lambda = c ./ f;
omegat1 = tau1*omega;
omegat2 = tau2*omega;
alpham = ((2*alphav1g*omegat1)./(lambda.*(1+(omegat1.^2)))) + ...
  ((2*alphav2g*omegat2)./(lambda.*(1+(omegat2.^2))));

alphacm = alpham + alphac;
