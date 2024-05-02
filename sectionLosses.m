function [G, Zc] = sectionLosses( r1, r2, L, f, T, lossType, humidity, P0 )
% SECTIONLOSSES: Compute the lossy complex wave propagation variable and
% characteristic impedance for a waveguide section.
%
% [G, ZC] = SECTIONLOSSES( R1, R2, L, F, T, LOSSYPE, HUMIDITY, P0 )
% computes the wave propagation variable GAMMA and the complex
% characteristic impedance ZC for given frequency values F in a
% cylindrical, conical or tonehole waveguide section of input radius R1,
% output radius R2 and length L. For cylindrical sections, R1 = R2. If the
% temperature T is not specified, a default value of 20 C is assumed. If
% either the relative HUMIDITY (as percentage) or atmospheric pressure (P0,
% in Pascals) is not specified, values of 40% and 101325 Pascals will be
% used by default. If parameter LOSSTYPE = 0, losses will be ignored.
% Otherwise, classical and molecular air losses are included, as well as
% various boundary layer loss approximations, with increasing accuracy as
% given by LOSSTYPE (default = 1):
%
%   1 - "Standard" TMM boundary layer losses from [1], [2], [3].
%       Zeroth/First order approximations.
%   2 - Second/Fourth order approximations (depending on bore radius)
%       from [4] and [5].
%   3 - Numerical computation of Bessel functions of thie first kind,
%       from [4].
%
% For LOSSTYPE = 1 and LOSSTYPE = 2, the type of loss (wide or narrow pipe)
% will change based on the maximum Stokes number, rv = R*sqrt(k/lv). If
% max(rv) < 1 then the narrow pipe approximation will be used.
%
% by Champ Darabundit and Gary Scavone, McGill University, 2023-2024.
%
% References:
%
% 1. A. Lefebvre, G.P. Scavone, and J. Kergomard (2013), "External tonehole
%       interactions in woodwind instruments." Acta Acustica united with
%       Acustica, vol. 99, pp. 975-985.
%
% 2. J. Abel, T. Smyth, and J.O. Smith (2003), "A simple, accurate wall
%       loss filter for acoustic tubes." in proc. of the 6th Int. Conf.
%       Digital Audio Effects (DAFx-03), London, UK.
%
% 3. A.H. Benade (1968), "On the propagation of sound waves in a
%       cylindrical conduit." Journal of the Acoustical Society of America,
%       vol. 44, pp. 616-623.
%
% 4. C. Zwikker and C.W. Kosten (1949), "Sound absorbing materials."
%       Elsevier Publishing Company.
%
% 5. A. Chaigne and J. Kergomard (2016), "Acoustics of music instruments."
%       Springer-Verlag.

if nargin < 4 || nargin > 8
  error( 'Invalid number of arguments.');
end
if ~isvector(f)
  error( 'f should be a 1D vector of frequencies in Hertz.' );
end
if ~exist( 'T', 'var')
  T = 20;
end
if ~exist( 'lossType', 'var')
  lossType = 1;
end
if ~exist( 'humidity', 'var')
  humidity = 40;  % room humidity as percentage
end
if ~exist( 'P0', 'var')
  P0 = 101325;    % atmospheric pressure at sea level in Pascals
end
if lossType > 3 || mod(lossType, 1) ~= 0
  error('lossType must be an integer between 0 and 3')
end

[c, rho, gamma, lv, Pr] = thermoConstants( T );
k = 2 * pi * f / c;
Zc = rho * c / ( pi * r1 * r2 );

if r1 == r2 % cylindrical or tonehole section
  req = r1;
else % conic section
  x1 = r1 * L / (r2 - r1);
  req = L * r1 / x1 / log(1 + L/x1);
end

% Case 0: No losses
if lossType == 0
  G = 1j*k;
  return;

% Case 1: TMM approximation from [1], [2], [3]
elseif lossType == 1
  rv = req*sqrt(abs(k/lv));
  if max(rv) < 1  % narrow pipe
    % G = k.*2.*sqrt(gamma)./rv + 1j.*k./rv.*2.*sqrt(gamma);
    G = k.*2.*sqrt(gamma)./rv.*(1 + 1j);
  else  % wide pipe
    cst = sqrt(lv/2).*(1 + (gamma - 1)/sqrt(Pr));
    G = 1j*k + (1 + 1j) * cst .* sqrt(k)/req;
  end

% Case 2: 2nd/4th order expansion approximations from [4], [5]
elseif lossType == 2
  lt = lv/Pr;
  rv = req*sqrt(abs(k/lv));
  rt = req*sqrt(abs(k/lt));
  S = pi*r1*r2;
  w = k.*c;
  if max(rv) < 1  % narrow pipe, this condition could be changed ...
    Zv = (rho*c/S).*(8*lv)./req^2.*(1 + (1/6).*1j.*rv.^2);
    Yt = 1j.*w.*gamma.*(S./(rho*c^2)).*(1 - (gamma - 1)./(8*gamma).*1j.*rt.^2);
  else  % wide pipe
    Zv = (1j*w*rho./S).*(1 + (2/sqrt(2))*(1 - 1j)./rv - 3j./rv.^2);
    Yt = (1j*w*S)./(rho*c^2).*(1 + (gamma - 1).*((2/sqrt(2))*(1 - 1j)./rt + 1j./rt.^2));
  end
  G = sqrt(Zv.*Yt);
  Zc = sqrt(Zv./Yt);

% Case 3: Numerical computation of Bessel function from [4], [5]
elseif lossType == 3
  lt = lv/Pr;
  kv = ((1 - 1j)/sqrt(2))*sqrt(k/lv);
  kt = ((1 - 1j)/sqrt(2))*sqrt(k/lt);
  S = pi*r1*r2;
  w = k.*c;
  % Compute Bessel functions with scaling in case of overflow
  J1v = besselj(1, kv.*req, 1); % J1v = J1v./ exp(-abs(imag(kv.*req)));
  J0v = besselj(0, kv.*req, 1); % J0v = J0v./ exp(-abs(imag(kv.*req)));
  J1t = besselj(1, kt.*req, 1); % J1t = J1t./ exp(-abs(imag(kt.*req)));
  J0t = besselj(0, kt.*req, 1); % J0t = J0t./ exp(-abs(imag(kt.*req)));
  % Compute immittances
  Zv = (1j.*w.*rho/S).*(1 - 2./(kv.*req).*(J1v./J0v)).^(-1);
  Yt = (1j.*w.*S)/(rho.*c^2).*(1 + (gamma - 1).*(2./(kt.*req)).*(J1t./J0t));
  G = sqrt(Zv.*Yt);
  Zc = sqrt(Zv./Yt);
end

% Include classical and molecular air losses
G = G + airLosses( f, T, humidity, P0 );

end