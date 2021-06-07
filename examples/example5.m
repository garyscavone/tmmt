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

% Get first geometry data
fingering = 0;
[boreData, holeData] = sevenCylinders( fingering );
if isempty( boreData )
  return;
end

% Do first TMM calculations and plot
plotTypes = [1 5];
Zin = tmm( boreData, holeData, rho, c, k, alpha, endType ); % cylinders
rzplot( f, Zin, plotTypes, true, false, [], 'r-');

% Get second geometry data
fingering = 0;
[boreData, holeData] = sevenSegments( fingering );
if isempty( boreData )
  return;
end

% Do second TMM calculations and plot
Zin = tmm( boreData, holeData, rho, c, k, alpha, endType ); % cones & cylinders
rzplot( f, Zin, plotTypes, true, true, [], 'b-'); % plot with initial hold on
legend('Cylinders only', 'Cylinders / Cones');
subplot(numel(plotTypes), 1, 1)
title('Input Impedance of 7 segment structure (lossy, Z_L unflanged)')
ylim([-40 40])