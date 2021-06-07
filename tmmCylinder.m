function [A, B, C, D] = tmmCylinder( k, L, r, Zc, alpha )
% TMMCYLINDER:  Compute the transfer-matrix coefficients for a cylindrical section.
%
% [A B C D Z0] = TMMCYLINDER(K, L, R, ZC, ALPHA) returns the coefficients of
% the transfer matrix describing a 1D cylindrical waveguide section. K is a
% vector of wave numbers (frequencies) at which the coefficients are
% computed, L is the section length, R is the section radius, and ZC is the
% characteristic impedance of the section. The optional parameter ALPHA is
% used to apply thermo-viscous losses. The returned values are column
% vectors with as many rows as elements in K.
%
% by Gary P. Scavone, McGill University, 2013-2020.

if ~isvector(k)
  error( 'k should be a 1D vector.' );
end

Gamma = 1j*k;
if ( nargin == 5 )
  if sum(size(k) ~= size(alpha))
    error( 'k and alpha should be the same size.' );
  end
  % Recompute Gamma for given radius
  Gamma = Gamma + (1+1j) .* alpha / r;
end

sinhL = sinh( L * Gamma );
coshL = cosh( L * Gamma );
A = coshL;
B = Zc * sinhL;
C = sinhL / Zc;
D = A;

