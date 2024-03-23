% A matlab script used to compute the theoretical input impedance or
% reflectance of a air column structure (as defined in a separate geometry
% file) using the TMM approach.
%
% by Gary P. Scavone, McGill University, 2021-2024.

clear; clf;
lossy = 1;   % 0 = lossless, 1 = traditional losses, 2 = Zwikker-Kosten; 3 = Bessel function
endType = 1; % 0 = closed, 1 = unflanged, 2 = flanged, 3 = ideally open

% Evaluation frequencies
fmax = 6000;          % maximum evaluation frequency (Hz)
N = fmax;             % number of frequencies for evaluation (even)
finc = fmax / (N-1);
f = eps:finc:fmax;
T = 20;   % temperature (C)

% Include path to needed scripts
addpath( '../', '../geometries/' );

% Get geometry data
fingering = 6;  % use a fingering with more open holes to see differences
[boreData, holeData] = keefeFlute( fingering );
if isempty( boreData )
  return;
end

% Do TMM calculations and plot
Zin = tmm( boreData, holeData, endType, f, lossy, T ); % tmm
rzplot( f, Zin, [1 11], true, false, [], 'r-');

Zin = tmmi( boreData, holeData, endType, f, lossy, T ); % tmmi
rzplot( f, Zin, [1 11], true, true, [], 'b-'); % plot with initial hold on
xlim([0 10]);
subplot(2, 1, 1)
ylim([-40 40])
title('Input Impedance / reflection function for Keefe flute')
legend('TMM (No Interactions)', 'TMMI (Interactions)');
