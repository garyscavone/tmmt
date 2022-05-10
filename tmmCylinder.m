function [A, B, C, D] = tmmCylinder( k, L, r, Zc, cst, alphacm )
% TMMCYLINDER:  Compute the transfer-matrix coefficients for a cylindrical section.
%
% [A B C D] = TMMCYLINDER(K, L, R, ZC, CST, ALPHACM) returns the
% coefficients of the transfer matrix describing a 1D cylindrical waveguide
% section. K is a 1D vector of wave numbers (frequencies) at which the
% coefficients are computed, L is the section length, R is the section
% radius, and ZC is the real characteristic impedance of the section. The
% optional parameter CST, a wall loss constant that depends on the
% properties of air, is used to apply thermo-viscous losses. ALPHACM is an
% optional 1D vector (the same size as K) of attenuation values
% corresponding to molecular and classical losses in air. The returned
% values are vectors of the same dimension as K.
%
% by Gary P. Scavone, McGill University, 2013-2022.

if ~isvector(k)
  error( 'k should be a 1D vector.' );
end

Gamma = 1j*k;
if exist( 'cst', 'var') && cst > 0
  % Include wall losses
  Gamma = Gamma + (1+1j) * cst .* sqrt(k) / r;
end

if exist( 'alphacm', 'var')
  if size(k) ~= size(alphacm)
    error( 'Incompatible sizes of k and alphacm vectors.' );
  end
  Gamma = Gamma + alphacm;
end

sinhL = sinh( L * Gamma );
coshL = cosh( L * Gamma );
A = coshL;
B = Zc * sinhL;
C = sinhL / Zc;
D = A;

