function [A, B, C, D] = tmmCylCone( r1, r2, L, Gamma, Zc )
% TMMCYLCONE:  Compute the transfer-matrix coefficients for a cylindrical
% or conical section.
%
% [A B C D] = TMMCYLCONE( R1, R2, L, GAMMA, ZC ) returns the coefficients
% of the transfer matrix describing a 1D cylindrical or conical waveguide
% section. R1 is the section input radius, R2 is the section output radius,
% L is the section length, GAMMA is a 1D vector of wave propagation values
% (wave numbers, with or without losses) at which the coefficients are
% computed, and ZC is the characteristic impedance of the section. For a
% cylindircal section, R1 and R2 should be the same values. The returned
% values are vectors of the same dimension as GAMMA.
%
% Initially by Gary P. Scavone, McGill University, 2013-2024, updates
% provided by Champ Darabundit and Miranda Jackson, 2023-2024.

if nargin ~= 5
  error( 'Invalid number of arguments.');
end
if ~isvector(Gamma)
  error( 'Gamma should be a 1D vector.' );
end

y1 = (r2 - r1)/(r1 * L); 
y2 = (r2 - r1)/(r2 * L);

sinhL = sinh( L * Gamma );
coshL = cosh( L * Gamma );
A = r2*coshL/r1 - y1*sinhL./Gamma;
B = Zc .* sinhL;
C = ((1 - y1*y2./(Gamma.^2)).*sinhL + ...
    (y1 - y2)*coshL./Gamma)./Zc;
D = r1*coshL/r2 + y2*sinhL./Gamma;
