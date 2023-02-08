% A matlab script used to compute the theoretical input impedance or
% reflectance of a air column structure (as defined in a separate geometry
% file) using the TMM approach.
%
% by Gary P. Scavone, McGill University, 2021-2022.

%clear; clf;
lossy = true; % turn on/off losses
endType = 1;   % 0 = closed, 1 = unflanged, 2 = flanged, 3 = ideally open

% Evaluation frequencies
fmax = 10000;         % maximum evaluation frequency (Hz)
N = fmax;             % number of frequencies for evaluation (even)
finc = fmax / (N-1);
f = eps:finc:fmax;
T = 20;   % temperature (C)

% Include path to needed scripts
addpath( '../', '../geometries/' );

% Get  geometry data
figure(1)
fingering = 3;
drawBore 'keefeFlute';
[boreData, holeData] = keefeFlute( fingering );
if isempty( boreData )
  return;
end

% Do TMM calculations and plot
figure(2)
plotTypes = [1 11];
Zin = tmm( boreData, holeData, f, T, lossy, endType );
rzplot( f, Zin, plotTypes, true, false, [], 'b-', true); % with time-domain smoothing
xlim([0 10]);
subplot(numel(plotTypes), 1, 1)
title('Input Impedance (top) and reflection function (bottom) for Keefe flute.')
ylim([-40 40])