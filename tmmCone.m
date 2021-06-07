function [A, B, C, D] = tmmCone( k, L, r1, r2, Zc, alpha )
% TMMCONE:  Compute the transfer-matrix coefficients for a conical section.
%
% [A B C D] = TMMCONE(K, L, R1, R2, ZC, ALPHA) returns the coefficients
% of the transfer matrix describing a 1D conical waveguide section. K
% is a vector of wave numbers (frequencies) at which the coefficients are
% computed, L is the section length, R1 is the section input radius, R2 is
% the section output radius, and ZC is the real characteristic impedance at
% the input. The optional parameter ALPHA is used to apply thermo-viscous
% losses. The returned values are column vectors with as many rows as
% elements in K.
%
% by Gary P. Scavone, McGill University, 2013-2021.

if ~isvector(k)
  error( 'k should be a 1D vector.' );
end

% Conical dimensions
x1 = r1 * L / (r2 - r1);
x2 = x1 + L;
req = L * r1 / x1 / log(1 + L/x1);

Gamma = 1j*k;
if ( nargin == 6 )
  if sum(size(k) ~= size(alpha))
    error( 'k and alpha should be the same size.' );
  end
  % Recompute Gamma for input radius
  Gamma = Gamma + (1+1j) * alpha / req;
end

% Rescale characteristic impedance
Zc = r1 * Zc / r2;

kc = -1j * Gamma;
sinL = sin(L*kc);
cosL = cos(L*kc);
A = (r2/r1)*cosL - sinL./(k*x1);
B = 1j*Zc*sinL;
C = (1j*(1 + 1./(x1*x2*k.^2)).*sinL + ...
    (1/x1 - 1/x2)*cosL/1j./k)/Zc;
D = r1*cosL/r2 + sinL/x2./k;
