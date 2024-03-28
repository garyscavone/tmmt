% A matlab script used to compute the theoretical input impedance or
% reflectance of a air column structure (as defined in a separate geometry
% file) with different unflanged approximations using the TMM approach.
%
% by Champ Darabundit, McGill University, 2024.

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
figure(1)
drawBore 'pipe'
fingering = 0;
[boreData, holeData] = pipe( fingering );
if isempty( boreData )
  return;
end

% Do TMM calculations for unflanged end condition
% Here we use the default unflanged approximation 'dalmont'
ZinUfDalmont = tmm( boreData, holeData, endType, f, lossy, T ); 

% Use the radiation function to access other approximations
rEnd = boreData(2, end);    % Grab radius at end from boreData
ZlCausse = radiation(rEnd, f, T, 'causse');
ZlLs = radiation(rEnd, f, T, 'unflanged'); % Levine & Schwinger

% Call tmm with vector of load impedance
ZinUfCausse = tmm( boreData, holeData, ZlCausse, f, lossy, T );
ZinUfLs = tmm( boreData, holeData, ZlLs, f, lossy, T );
% Plot result
figure(2)
plot( f, 20*log10(abs(ZinUfDalmont)), 'b-', ...
      f, 20*log10(abs(ZinUfCausse)), 'r--', ...
      f, 20*log10(abs(ZinUfLs)), 'g:', ...
      'LineWidth', 2 );
ylabel('20*log10(|Impedance|)')
xlabel('Frequency (Hz)')
ylim([-100 100])
grid
title('Input Impedance of 60 cm pipe (Z_L based on different unflanged approximations)')
legend('Unflanged Dalmont', 'Unflanged Causse', 'Unflanged Levine & Schwinger');
