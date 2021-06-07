function [W, X, Y, Z] = tonehole( k, delta, b, t, rPad, hPad, w, state, Zc, chimney, alpha, type )
% TONEHOLE:  Compute the transfer matrix coefficients for a tonehole.
%
% [W X Y Z] = TONEHOLE(K, DELTA, B, T, RPAD, HPAD, W, STATE, ZC, CHIMNEY,
% ALPHA, TYPE) returns the coefficients of the transfer matrix describing a
% tonehole section. K is a vector of wave numbers (frequencies) at which
% the coefficients are computed, DELTA is the ratio of hole to air column
% radii, B is the hole radius, T is the height, RPAD is the radius of a
% hanging pad over the hole, HPAD is the distance of the pad from the hole,
% W is the wall thickness of the hole, STATE is 0 or 1 to indicate whether
% the hole is closed or open, respectively,  ZC is the characteristic
% impedance of the hole, and CHIMNEY is the length the tonehole extends out
% from the bore. The optional parameters ALPHA and TYPE are used to apply
% thermo-viscous losses and specify a particular characterization (default
% is 'Dalmont2002'). The returned values are column vectors with as many
% rows as elements in K.
%
% by Gary P. Scavone, McGill University, April 2013.
% Based in part on functions from WIAT by Antoine Lefebvre.
%
% References:
%
% 1. Dalmont, J. et al., "Experimental Determination of the Equivalent
%      Circuit of an Open Side Hole: Linear and Non Linear Behaviour,"
%      Acta Acustica united with Acustica,  vol. 88, pp. 567-575, 2002.
%
% 2. van Walstijn, M., & Campbell, M. (2003), "Discrete-time modeling of
%      woodwind instrument bores using wave variables," J. Acoust. Soc.
%      Am., vol. 113(1), 575-585.
%
% 3. Nederveen, C.J., Jansen, J.K.M, and van Hassel, R.R. (1998).
%      "Corrections for woodwind tone-hole calculations.",
%      Acta Acustica united with Acustica, vol. 84, pp. 957-966.
%
% 4. Lefebvre, A. and Scavone, G. (2012), "Characterization of woodwind
%      instrument toneholes with the finite element method.", J. Acoust.
%      Soc. Am., vol 131(4), pp. 3153-3163.

[M, N] = size( k );
if ( M == 1 )
  k = k.';
elseif ( M > 1 && N > 1 )
  error( 'k should be a 1D vector.' );
end

%if ( nargin < 12 )
%  type = 'Dalmont2002';
%end

if ( nargin >= 11 )
  % Recompute Gamma for given radius
  Gamma = 1j*k + (1+1j) * alpha / b;
  k = Gamma / 1j;
end

% Inner length correction: Eq. (4) in [1]
ti = b * ( 0.82 - 1.4*delta^2 + 0.75*delta^2.7 );

% Matching volume length correction: Eq. (6) in [1]
tm = b * delta * (1. + 0.207*delta^3) / 8;

if state % open hole
  % Series length correction: Eq. 33) in [4]
  ta = (-0.35 + 0.06*tanh(2.7*t/b)) * b * delta^2;

  % Pad correction: Eq. (48) in [1]
  tp = 0;
  if rPad > 0
    tp = b / (3.5*(hPad/b)^0.8*(hPad/b+3*w/b)^(-0.4)+30*(hPad/rPad)^2.6);
  end

  % Open hole radiation impedance (low-frequency approximation)
  ZroZc = 0.25*(k*b).^2 + 1j*k*(0.6113*b + tp );

  if chimney < 0.002
    %ZroZc = ZroZc + 0.25*(k*b).^2 + 1j*k*(0.8215-0.6113)*b;
  end

  % Radiation length correction
  tr = atan(-1j*ZroZc)./k;

  % Shunt impedance
  Zs = 1j * Zc * (k.*ti + tan( k.*(t+tm+tr) ));

else % close hole
  % Series length correction: Eq. (34) in [4]
  ta = (-0.12 - 0.17*tanh(2.4*t/b)) * b * delta^2;

  % Shunt impedance
  Zs = 1j * Zc * (k*ti - cot( k * (t + tm) ));
end

% Series impedance
Za = 1j * Zc * delta^2 * k * ta;

ZaoZs = Za ./ Zs;
W = 1 + 0.5 * ZaoZs;
X = Za .* (1 + 0.25 * ZaoZs);
Y = 1 ./ Zs;
Z = W;

if ( nargin < 12 ), return; end

if strcmp( type, 'KeefeMatrix' )
  W = 1; X = Za; Z = 1;
end

