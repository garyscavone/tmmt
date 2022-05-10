% A matlab script used to compute the theoretical input impedance or
% reflectance of a air column structure (as defined in a separate geometry
% file) using the TMM approach.
%
% by Gary P. Scavone, McGill University, 2021-2022.

clear; clf;
lossy = false; % turn on/off losses
endType = 3;   % 0 = closed, 1 = unflanged, 2 = flanged, 3 = ideally open

% Evaluation frequencies
fmax = 6000;          % maximum evaluation frequency (Hz)
N = fmax;             % number of frequencies for evaluation (even)
finc = fmax / (N-1);
f = eps:finc:fmax;
T = 20;   % temperature (C)

% Include path to needed scripts
addpath( '../', '../geometries/' );

% Get geometry data
figure(1)
drawBore 'pipe'
fingering = 0;
[boreData, holeData] = pipe( fingering );
if isempty( boreData )
  return;
end

% Do TMM calculations
ZinOpen = tmm( boreData, holeData, f, T, lossy, endType ); % ideally open end
ZinClosed = tmm( boreData, holeData, f, T, lossy, 0 );     % ideally closed end

% Plot result
figure(2)
plot( f, 20*log10(abs(ZinOpen)), 'b-', f, 20*log10(abs(ZinClosed)), 'r-', 'LineWidth', 2 );
ylabel('20*log10(|Impedance|)')
xlabel('Frequency (Hz)')
ylim([-100 100])
grid
title('Input Impedance of 60 cm pipe (Z_L = 0 and Z_L = \infty)')
legend('Open End', 'Closed End');
