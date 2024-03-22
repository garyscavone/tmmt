function Zin = tmm( boreData, holeData, endType, f, lossType, T )
% TMM: Compute the normalized input impedance of a system using the
%      transfer matrix method.
%
% ZIN = TMM( BOREDATA, HOLEDATA, ENDTYPE, F, LOSSTYPE, T ) returns the
% input impedance of a system defined by BOREDATA and HOLEDATA, normalized
% by the characteristic impedance at the input, at frequencies specified in
% the 1D vector F, given an optional air temperature T in degrees Celsius
% (default = 20 C). The parameter ENDTYPE specifies the bore end condition
% [0 = rigidly closed; 1 = unflanged open; 2 = flanged open; 3 = ideally
% open (Zl = 0)]. The optional parameter LOSSTYPE specifies how losses are
% approximated [0 = no losses; 1 = lowest order losses (previous tmm
% method, default); 2 = Zwikker-Kosten; 3 = full Bessel function
% computations].
%
% Initially by Gary P. Scavone, McGill University, 2013-2024, updates
% provided by Champ Darabundit, 2023.

if nargin < 4 || nargin > 6
  error( 'Invalid number of arguments.');
end
if ~isvector(f)
  error( 'f should be a 1D vector of frequencies in Hertz.' );
end
if ~exist( 'T', 'var')
  T = 20;
end
if ~exist( 'lossType', 'var')
  lossType = 1;
end

% Bore dimensions
idx = find(diff(boreData(1,:)) == 0);
boreData(1, idx+1) = boreData(1, idx+1) + eps; % avoid double values
x = sort( [boreData(1,:) holeData(1,:)] ); % segment positions along x-axis
L = diff( x );                             % lengths of segments

isHole = zeros(size(x));                   % is x value at a tonehole?
for n = 1:length(x)
  isHole(n) = 1 - isempty(find(x(n)==holeData(1,:), 1));
end

% Interpolate bore radii at x values
ra = interp1(boreData(1,:), boreData(2,:), x, 'linear');

% Tonehole dimensions and states
rb = holeData(2,:).';            % tonehole radii
t = holeData(3,:).';             % tonehole heights
chimney = holeData(4,:);         % tonehole chimney height
states = holeData(5,:);          % tonehole states
[n, m] = size(holeData);
padr = zeros(1, m);
padt = zeros(1, m);
holew = zeros(1, m);
if n > 6, padr = holeData(7,:); end % tonehole pad radii
if n > 7, padt = holeData(8,:); end % tonehole pad heights
if n > 8, holew = holeData(9,:); end % tonehole wall thickness

% Work our way back from the load impedance at the end.
switch endType
  case 1
    Zl = radiation( ra(end), f, T, 'dalmont' ); % L&S unflanged approximation
  case 2
    Zl = radiation( ra(end), f, T, 'flanged' ); % load impedance at end
  case 3
    Zl = 0;
  otherwise
    Zl = Inf;
end

nHole = sum(isHole);
for n = length(L):-1:1
  if L(n) > eps
    [Gamma, Zc] = sectionLosses( ra(n), ra(n+1), L(n), f, T, lossType );
    [A, B, C, D] = tmmCylCone( ra(n), ra(n+1), L(n), Gamma, Zc );
    if Zl == 0
      Zl = B ./ D;
    else
      Zl = (A + B./Zl) ./ (C + D./Zl);
    end
  end
  
  if isHole(n)
    [Gamma, ~] = sectionLosses( rb(nHole), rb(nHole), 0, f, T, lossType );
    [A, B, C, D] = tmmTonehole( rb(nHole)/ra(n), rb(nHole), t(nHole), ...
      states(nHole), Gamma, '', T, chimney(nHole), padr(nHole), ...
      padt(nHole), holew(nHole) );
    nHole = nHole - 1;
    Zl = (A + B./Zl) ./ (C + D./Zl);
  end
end

Zin = Zl ./ Zc;
