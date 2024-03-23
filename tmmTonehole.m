function [W, X, Y, Z] = tmmTonehole( delta, b, height, state, Gamma, type, T, chimney, rPad, hPad, w )
% TMMTONEHOLE:  Compute the transfer matrix coefficients for a tonehole.
%
% [W X Y Z] = TMMTONEHOLE(DELTA, B, HEIGHT, STATE, GAMMA, TYPE, T, CHIMNEY,
% RPAD, HPAD) returns the coefficients of the transfer matrix describing a
% tonehole section. DELTA is the ratio of hole to air column radii, B is
% the hole radius, HEIGHT is the tonehole height, STATE is 0 or 1 to
% indicate whether the hole is closed or open, and GAMMA is a 1D vector of
% wave propagation values (wave numbers, with or without losses) at which
% the coefficients are computed. The remaining parameters are optional:
% TYPE is either 'Keefe1990', 'Dalmont2002', or 'Lefebvre2012' (default) to
% specify the model, T is the air temperature in degrees Celsius (default =
% 20 C), CHIMNEY is the length the tonehole extends out from the bore, RPAD
% is the radius of a hanging pad over the hole (use 0 for no pad), HPAD is
% the distance of the pad from the hole and W is the wall thickness of the
% tonehole (default values = 0 if not otherwise specified). The returned
% values are vectors of the same dimension as GAMMA.
%
% by Gary P. Scavone, McGill University, 2013-2024.
% Based in part on functions from WIAT by Antoine Lefebvre.
%
% References:
%
% 1. Keefe, D. H. (1990), "Woodwind air column models." J. Acoust. Soc.
%      Am., vol. 88(1), pp. 35–51.
%
% 2. Nederveen, C.J., Jansen, J.K.M, and van Hassel, R.R. (1998),
%      "Corrections for woodwind tone-hole calculations." Acta Acustica
%      united with Acustica, vol. 84, pp. 957-966.
%
% 3. Dubos, V., Kergomard, J., Keefe, D., Khettabi, A., Dalmont, J.P., 
%      Nederveen, C.J. (1999), "Theory of sound propagation in a duct with
%      a branched tube using modal decomposition." Acta Acustica united
%      with Acustica, vol. 85, 153–169.
%
% 4. Dalmont, J.P. and Nederveen, C.J. (2001), "Radiation impedance of
%      tubes with different flanges: numerical and experimental
%      investigations." J. Sound and Vib., vol. 244, pp. 505-534.
%
% 5. Dalmont, J.P., Nederveen, C.J., Dubos, V., Ollivier, S., Meserette,
%      V., Sligte, E. (2002), "Experimental determination of the equivalent
%      circuit of an open side hole: Linear and non linear behaviour." Acta
%      Acustica united with Acustica, vol. 88, pp. 567-575.
%
% 6. Lefebvre, A. and Scavone, G. (2012), "Characterization of woodwind
%      instrument toneholes with the finite element method.", J. Acoust.
%      Soc. Am., vol 131(4), pp. 3153-3163.

if nargin < 5 || nargin > 11
  error( 'Incorrect number of parameters.' );
end
if ~isvector(Gamma)
  error( 'Gamma should be a 1D vector.' );
end
if ~exist( 'T', 'var')
  T = 20;
end

[c, rho, ~, lv, ~] = thermoConstants( T );
k = -1j * Gamma;
Zc = rho * c / ( pi * b * b );

if nargin < 8 || isempty(type), type = 'Lefebvre2012'; end
if nargin < 9, chimney = 0; end
if nargin < 10, rPad = 0; end
if nargin < 10, hPad = 0; end
if nargin < 11, w = 0; end

if strcmp( type, 'Keefe1990' )
  tgh = height + b*delta*(1 + 0.172*delta^2)/8; % tonehole geometric height

  if state % open hole
    dv = sqrt(lv ./ k);
    rc = 0.0005; % tonehole radius of curvature (m) ... estimate

    if rPad > 0 % effective length with pad
      tmp = 0.61*((rPad/b)^0.18)*((b/hPad)^0.39);
      te = ((1./k).*tan(tgh*k) + b*(tmp + (pi/4)*(1 - 0.74*delta^2))) ./ ...
        (1 - tmp*b*k.*tan(tgh*k));
    else        % effective length without pad
      te = (tan(tgh*k)./k + b*(1.4 - 0.58*delta^2)) ./ ...
        (1 - 0.61*b*k.*tan(tgh*k));
    end
    
    % Specific resistance
    xie = 0.25*(b*k).^2 + tgh*imag(k) + 0.25*log(2*b/rc)*dv.*k;

    % Open series equivalent length
    tao = (0.47*b*delta^4)/(tanh(1.84*tgh/b) + 0.62*delta^2 + 0.64*delta);

    % Open tonehole shunt & series impedances
    Zs = Zc .* (1j*k.*te + xie);
    Za = -1j * Zc .* k * tao;

  else % closed hole
    % Closed series equivalent length
    tac = (0.47*b.*delta.^4)/(coth(1.84*tgh/b) + 0.62*delta^2 + 0.64*delta);
    
    % Closed tonehole shunt & series impedances
    Zs = -1j*Zc .* (cot(tgh*k) + tgh*(0.25*(b./tgh).^2 ...
      + 0.58*delta^2 - 0.25*pi*b / tgh)*k);
    Za = -1j * Zc .* k * tac;
  end
  W = 1; X = Za; Y = 1 ./ Zs; Z = 1;
  return;

elseif strcmp( type, 'Dalmont2002' )
  % Inner length correction: Eq. (40) in [2]
  ti = b * ( 0.82 - 1.4*delta^2 + 0.75*delta^2.7 );

  if state % open hole series length correction from [2]
    ta = -0.28 * b * delta^2;
  else % closed hole series length correction from [3]
    ta = -b * delta^2 / (1.78*coth(1.84*height/b) + 0.94 + 0.54*delta + ...
      0.285*delta^2);
  end
  
elseif strcmp( type, 'Lefebvre2012' )
  % Frequency-dependent inner length correction: Eqs. (31-32) in [6]
  ti = b * ( 0.822 - 0.095*delta - 1.566*delta^2 + 2.138*delta^3 ...
    - 1.64*delta^4 + 0.502*delta^5 );
  H = 1 - 4.56*delta + 6.55*delta^2;
  ka = k * b / delta;
  I = 0.17*ka + 0.92*ka.^2 + 0.16*ka.^3 - 0.29*ka.^4;
  ti = ti * (1 + H*I);

  if state % open hole series length correction: Eq. (33) in [6]
    ta = (-0.35 + 0.06*tanh(2.7*height/b)) * b * delta^2;
  else % closed hole series length correction: Eq. (34) in [6]
    ta = (-0.12 - 0.17*tanh(2.4*height/b)) * b * delta^2;
  end
  
else
  error( 'Unknown model type parameter.' );
end

% Matching volume length correction from [2]
tm = b * delta * (1. + 0.207*delta^3) / 8;

if state % open hole

  % Pad correction: Eq. (48) in [4]
  tp = 0;
  if rPad > 0
    tp = b / (3.5*(hPad/b)^0.8*(hPad/b+3*w/b)^(-0.4)+30*(hPad/rPad)^2.6);
  end

  % Open hole radiation impedance (unflanged low-frequency approximation)
  ZroZc = 0.25*(k*b).^2 + 1j*k*(0.6113*b + tp );

%   if chimney < 0.002
%     ZroZc = ZroZc + 0.25*(k*b).^2 + 1j*k*(0.8215-0.6113)*b;
%   end

  % Radiation length correction
  tr = atan(-1j*ZroZc)./k;

  % Shunt impedance
  Zs = 1j * Zc .* (k.*ti + tan( k.*(height+tm+tr) ));

else % closed hole

  % Shunt impedance
  Zs = 1j * Zc .* (k.*ti - cot( k * (height + tm) ));

end

% Series impedance
Za = 1j * Zc .* delta^2 * k * ta;

ZaoZs = Za ./ Zs;
W = 1 + 0.5 * ZaoZs;
X = Za .* (1 + 0.25 * ZaoZs);
Y = 1 ./ Zs;
Z = W;
