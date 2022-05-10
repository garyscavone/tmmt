% A matlab script used to compute the theoretical input impedance or
% reflectance of a air column structure (as defined in a separate geometry
% file) using the TMM approach.
%
% by Gary P. Scavone, McGill University, 2021-2022.

clear; clf;
lossy = true; % turn on/off losses
endType = 1;   % 0 = closed, 1 = unflanged, 2 = flanged, 3 = ideally open

% Evaluation frequencies
fmax = 6000;          % maximum evaluation frequency (Hz)
N = fmax;             % number of frequencies for evaluation (even)
finc = fmax / (N-1);
f = eps:finc:fmax;
T = 20;   % temperature (C)

% Include path to needed scripts
addpath( '../', '../geometries/' );

% Get first geometry data
fingering = 0;
[boreData, holeData] = sevenCylinders( fingering );
if isempty( boreData )
  return;
end

% Do first TMM calculations and plot
figure(1)
plotTypes = [1 6];
Zin = tmm( boreData, holeData, f, T, lossy, endType ); % cylinders
rzplot( f, Zin, plotTypes, true, true, [], 'r-');

% Get second geometry data
fingering = 0;
figure(2)
drawBore 'sevenSegments';
[boreData, holeData] = sevenSegments( fingering );
if isempty( boreData )
  return;
end

% Do second TMM calculations and plot
figure(1)
Zin = tmm( boreData, holeData, f, T, lossy, endType ); % cones & cylinders
rzplot( f, Zin, plotTypes, true, false, [], 'b-'); % plot with initial hold on
legend('Cylinders only', 'Cylinders / Cones');
subplot(numel(plotTypes), 1, 1)
title('Input Impedance of 7 segment structure (lossy, Z_L unflanged)')
ylim([-40 40])