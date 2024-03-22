function alpha = airLosses( f, T, humidity, P0 )
%  AIRLOSSES: Determine classical and molecular propagation losses at the
%  specified frequencies for given values of temperature, humidity and
%  atmospheric pressure.
%
% ALPHA = AIRLOSSES( F, T, HUMIDITY, P0 ) returns a vector ALPHA (of size
% equal to F) of wave propagation loss factors due to classical and
% molecular effects. If the temperature T is not specified, a default value
% of 20 C is assumed. If either the relative HUMIDITY (as percentage) or
% atmospheric pressure (P0, in Pascals) is not specified, values of 40% and
% 101325 Pascals will be used by default.
%
% References:
%
% 1. M. van Walstijn, M. Campbell, J. Kemp, D. Sharp, "Wideband
%      Measurements of the Acoustic Impedance of Tubular Objects," Acta
%      Acustica united with Acustica, Vol. 91, pp. 590-604, 2005.
%
% by Gary P. Scavone, McGill University, 2013-2024.

if nargin < 1 || nargin > 4
  error( 'Invalid number of arguments.');
end
if ~isvector(f)
  error( 'f should be a 1D vector of frequencies in Hertz.' );
end
if ~exist( 'T', 'var')
  T = 20;
end
if ~exist( 'humidity', 'var')
  humidity = 40;  % room humidity as percentage
end
if ~exist( 'P0', 'var')
  P0 = 101325;    % atmospheric pressure at sea level in Pascals
end

[c, ~, gamma, lv, Pr] = thermoConstants( T );

omega = 2 * pi * f;
alphac = omega.^2 * lv * (4/3 + 0.6 + (gamma-1)/Pr) / (2*c^2); % classical attenuation

pv = 0.065773*T^3 + 0.1445*T^2 + 59.34*T + 560.54; % vapor pressure of water
cv1oR = (-0.000002778*T^3 + 0.0007857*T^2 + 0.08599*T + 3.883)*10^-3; % ratio of specific-heat to gas constant for O2
cv2oR = (-0.00000009259*T^3 + 0.00035596*T^2 + 0.02212*T + 0.5525)*10^-3; % ratio of specific-heat to gas constant for N2
alphavgm = pi*(gamma-1)^2/(2*gamma);
alphav1g = alphavgm * cv1oR;  % maximum absorption per wavelength for O2
alphav2g = alphavgm * cv2oR;  % maximum absorption per wavelength for N2

h = 0.01 * humidity * pv / P0;  % humidity factor
G = 4.41*10^6*h*(0.05 + 100*h)/(0.391 + 100*h);
T1 = 293.16;
TA = T + 273.15; % absolute temperature in Kelvin
F = 6.142*((T1/TA)^(1/3) - 1);
tau1 = 1/(2*pi*(24 + G)); % relaxation time for oxygen in air (21%)
tau2 = sqrt(TA/T1)/(2*pi*(9 + (3.5*10^4)*h*exp(-F))); % relaxation time for nitrogen in air (78%)

lambda = c ./ f;
omegat1 = tau1*omega;
omegat2 = tau2*omega;
alpham = ((2*alphav1g*omegat1)./(lambda.*(1+(omegat1.^2)))) + ...
  ((2*alphav2g*omegat2)./(lambda.*(1+(omegat2.^2))));

alpha = alpham + alphac;
