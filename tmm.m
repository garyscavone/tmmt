function Zin = tmm( boreData, holeData, f, T, lossy, endType )
% TMM: Compute the normalized input impedance of a system using the
%      transfer matrix method.
%
% ZIN = TMM( BOREDATA, HOLEDATA, F, T, LOSSY, ENDTYPE ) returns the input
% impedance of a system defined by BOREDATA and HOLEDATA, normalized by the
% characteristic impedance at the input, at frequencies specified in the 1D
% vector F, given an optional air temperature T in degrees Celsius (default
% = 20 C). If the optional parameter LOSSY = false, losses will be ignored
% (default = true). The optional parameter ENDTYPE specifies the bore end
% condition [0 = rigidly closed; 1 = unflanged open (default); 2 = flanged
% open; 3 = ideally open (Zl = 0)].
%
% by Gary P. Scavone, McGill University, 2013-2022.

if nargin < 3 || nargin > 6
  error( 'Invalid number of arguments.');
end
if ~isvector(f)
  error( 'f should be a 1D vector of frequencies in Hertz.' );
end
if ~exist( 'T', 'var')
  T = 20;
end
if ~exist( 'lossy', 'var')
  lossy = true;
end
if ~exist( 'endType', 'var')
  endType = 1;
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

% Get physical variables and attenuation values
[c, rho, wallcst, alphacm] = physicalSettings( T, f );
if ~lossy
  wallcst = 0;
  alphacm = alphacm * 0;
end
k = 2 * pi * f / c;

% Compute characteristic impedances for each bore and tonehole radius
Zc = rho * c ./ (pi * ra.^2);
Zch = rho * c ./ (pi * rb.^2);

% Work our way back from the load impedance at the end.
switch endType
  case 1
    Zl = Zc(end)*radiation( k, ra(end), 'dalmont' ); % L&S unflanged approximation
  case 2
    Zl = Zc(end)*radiation( k, ra(end), 'flanged' ); % load impedance at end
  case 3
    Zl = 0;
  otherwise
    Zl = Inf;
end

nHole = sum(isHole);
for n = length(L):-1:1
  if L(n) > eps
    if ( ra(n) == ra(n+1) )
      [A, B, C, D] = tmmCylinder( k, L(n), ra(n), Zc(n), wallcst, alphacm );
    else
      [A, B, C, D] = tmmCone( k, L(n), ra(n), ra(n+1), Zc(n), wallcst, alphacm );
    end
    if Zl == 0
      Zl = B ./ D;
    else
      Zl = (A + B./Zl) ./ (C + D./Zl);
    end
  end
  
  if isHole(n)
    [A, B, C, D] = tmmTonehole( k, rb(nHole)/ra(n), rb(nHole), t(nHole), ...
      Zch(nHole), states(nHole), wallcst, alphacm, '', chimney(nHole), padr(nHole), ...
      padt(nHole), holew(nHole) );
    nHole = nHole - 1;
    Zl = (A + B./Zl) ./ (C + D./Zl);
  end
end

Zin = Zl / Zc(1);
