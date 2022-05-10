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

% Get geometry data
fingering = 7;
[boreData, holeData] = keefeFlute( fingering );
if isempty( boreData )
  return;
end

% Do TMM calculations and plot
Zin = tmm( boreData, holeData, f, T, lossy, endType ); % tmm
rzplot( f, Zin, 1, true, false, [], 'r-');

Zin = tmmi( boreData, holeData, f, T, lossy, endType ); % tmmi
rzplot( f, Zin, 1, true, true, [], 'b-'); % plot with initial hold on
title('Input Impedance for Keefe flute (all holes open)')
legend('TMM (No Interactions)', 'TMMI (Interactions)');
ylim([-40 40])