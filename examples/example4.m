% A matlab script used to compute the theoretical input impedance or
% reflectance of a air column structure (as defined in a separate geometry
% file) using the TMM approach.
%
% by Gary P. Scavone, McGill University, 2021.

clear; clf;
lossy = true; % turn on/off losses
endType = 1;   % 0 = closed, 1 = unflanged, 2 = flanged, 3 = ideally open

% Evaluation frequencies
fmax = 6000;          % maximum evaluation frequency (Hz)
N = fmax;             % number of frequencies for evaluation (even)
finc = fmax / (N-1);
f = eps:finc:fmax;
omega = 2*pi*f.';

% Include path to needed scripts
addpath( '../', '../geometries/' );

% Physical constants
T = 20;   % temperature (C)
[c, rho, CST] = physicalSettings( T );
k = omega / c;
if ~lossy
  CST = 0;
end
alpha = sqrt(k) * CST;             % loss factor, not including radius

% Get geometry data
fingering = 0;
[boreData, holeData] = pipe( fingering );
if isempty( boreData )
  return;
end

% Do TMM calculations and plot
Zin = tmm( boreData, holeData, rho, c, k, alpha, endType ); % unflanged
rzplot( f, Zin, 1, true, false, [], 'ro');

Zin = tmm( boreData, holeData, rho, c, k, alpha, 3 ); % ideally open
rzplot( f, Zin, 1, true, true, [], 'b-'); % plot with initial hold on
title('Input Impedance of 60 cm pipe (Z_L = 0 vs. unflanged)')
legend('Unflanged', 'Ideally Open');
ylim([-40 40])