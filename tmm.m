function Zin = tmm( boreData, holeData, rho, c, k, alpha, endType )
% TMM: Compute the normalized input impedance of a system using the
%      transfer matrix method.
%
% ZIN = TMM( BOREDATA, HOLEDATA, RHO, C, K, GAMMA, ALPHA, ENDTYPE ) returns
% the input impedance of a system, normalized by the characteristic
% impedance at the input, described by BOREDATA and HOLEDATA at frequencies
% specified by the wavenumber K, given values of air mass density RHO,
% speed of sound C, and loss factor ALPHA. The optional parameter ENDTYPE
% specifies the bore end condition [0 = rigidly closed; 1 = unflanged open
% (default); 2 = flanged open; 3 = ideally open (Zl = 0)].
%
% by Gary P. Scavone, McGill University, 2013-2021.

if nargin < 7
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
padw = zeros(1, m);
if ( n > 6 )
  padr = holeData(7,:);          % tonehole pad radii
  padt = holeData(8,:);          % tonehole pad heights
  padw = holeData(9,:);          % tonehole wall thickness
end

% Compute characteristic impedances for each bore and tonehole radius
Zc = rho * c ./ (pi * ra.^2);
Zch = rho * c ./ (pi * rb.^2);

% Work our way back from the load impedance at the end.
switch endType
  case 1
    Zl = Zc(end)*radiation( k, ra(end), 'UnflangedDalmont' ); % L&S unflanged approximation
  case 2
    Zl = Zc(end)*radiation( k, ra(end), 'FlangedNorris' ); % load impedance at end
  case 3
    Zl = 0;
  otherwise
    Zl = Inf;
end

nHole = sum(isHole);
for n = length(L):-1:1
  if L(n) > eps
    if ( ra(n) == ra(n+1) )
      [A, B, C, D] = tmmCylinder( k, L(n), ra(n), Zc(n), alpha );
    else
      [A, B, C, D] = tmmCone( k, L(n), ra(n), ra(n+1), Zc(n), alpha );
    end
    if Zl == 0
      Zl = (A.*Zl + B) ./ (C.*Zl + D);
    else
      Zl = (A + B./Zl) ./ (C + D./Zl);
    end
  end
  
  if isHole(n)
    [A, B, C, D] = tonehole( k, rb(nHole)/ra(n), rb(nHole), t(nHole), ...
      padr(nHole), padt(nHole), padw(nHole), states(nHole), Zch(nHole), ...
      chimney(nHole), alpha );
    nHole = nHole - 1;
    Zl = (A + B./Zl) ./ (C + D./Zl);
  end
end

Zin = Zl / Zc(1);
