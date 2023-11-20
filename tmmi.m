function Zin = tmmi( boreData, holeData, f, T, lossy, endType )
% TMMI: Compute the normalized input impedance of a system using the
%       transfer matrix method with external tonehole interactions.
%
% ZIN = TMMI( BOREDATA, HOLEDATA, F, T, LOSSY, ENDTYPE ) returns the input
% impedance of a system defined by BOREDATA and HOLEDATA, normalized by the
% characteristic impedance at the input, at frequencies specified in the 1D
% vector F, given an optional air temperature T in degrees Celsius (default
% = 20 C). If the optional parameter LOSSY = false, losses will be ignored
% (default = true). The optional parameter ENDTYPE specifies the bore end
% condition [0 = rigidly closed; 1 = unflanged open (default); 2 = flanged
% open; 3 = ideally open (Zl = 0)].
%
% by Gary P. Scavone, McGill University, 2013-2022.
%
% References:
%
% 1. Lefebvre, A., Scavone, G. and Kergomard, J. (2013), "External Tonehole
%      Interactions in Woodwind Instruments." Acta Acustica united with
%      Acustica, vol. 99, pp. 975-985.
%
% 2. Kergomard, J. (1989), "Tone hole external interactions in woodwind
%      musical instruments." Proceedings of the 1989 Congress on Acoustics,
%      Belgrade.

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

% Calculate and use mutual radiation impedances. If this is off, the
% calculations should be equivalent to the TMM.
doInteractions = true;
  
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
if ( n > 6 )
  padr = holeData(7,:);          % tonehole pad radii
  padt = holeData(8,:);          % tonehole pad heights
  holew = holeData(9,:);         % tonehole wall thickness
end

nOth = sum( states );   % number of open toneholes

% Check pipe end condition
nOpen = nOth;
if endType
  nOpen = nOth + 1;
end

if nOpen < 2  % Do TMM
  Zin = tmm( boreData, holeData, f, T, lossy, endType );
  return
end

% Get physical variables and attenuation values
[c, rho, gamma, lv, Pr, alphacm] = physicalSettings( T, f );
if ~lossy
  wallcst = 0;
  alphacm = alphacm * 0;
end
k = 2 * pi * f / c;

% Compute characteristic impedances for each bore and tonehole radius
Zc = rho * c ./ (pi * ra.^2);
Zch = rho * c ./ (pi * rb.^2);

% Determine indices of open holes in position and tonehole vectors
oidx = find(states);           % open tonehole indices (relative to all toneholes)
xidx = find(isHole);           % x-indices of holes
xidx = [xidx(oidx) length(x)]; % x-indices of open tone holes + end
ds = holeData(1, oidx);        % positions of open holes
ras = ra(isHole>0);            % radii at all toneholes

% Allocate various matrices
Ymunm1 = 0;
Ypnm1 = 0;
ZB = zeros( nOpen, nOpen, length(k) );
Y = zeros( nOpen, nOpen, length(k) );

for n = 1:nOth

  nHole = oidx(n);
  [~, B, C, ~] = tmmTonehole( k, rb(nHole)/ras(nHole), rb(nHole), t(nHole), ...
    Zch(nHole), states(nHole), wallcst, alphacm, '', chimney(nHole), ...
    padr(nHole), padt(nHole), holew(nHole) );
      
  % Diagonals of (Z+B) matrix are open hole shunt impedances.
  ZB(n, n, :) = 1 ./ C;
  
  % Compute other Z terms (mutual radiation impedances)
  if doInteractions
    for m = 1:nOth
      if m == n, continue; end
      dnm = abs( ds(n) - ds(m) );
      %ZB(n, m, :) = 1j*rho*c*k.*exp(-1j*k*dnm)/(4*pi*dnm);
      ZB(n, m, :) = 1j*rho*c*k.*exp(-1j*k*dnm)/(2*pi*dnm); % extra factor of 2
    end
  end
    
  % Compute Y elements: Start with 1/2 of series length correction for
  % current open hole.
  MA = 1; MB = B/2; MC = 0; MD = 1;

  % Cascade all sections between current open hole and the next open hole.
  for m = xidx(n):xidx(n+1)-1
    if isHole(m)
      if states(nHole) == 0 % closed
        [A, B, C, D] = tmmTonehole( k, rb(nHole)/ras(nHole), rb(nHole), ...
          t(nHole), Zch(nHole), states(nHole), wallcst, alphacm, '', ...
          chimney(nHole), padr(nHole), padt(nHole), holew(nHole) );
        MAT = MA.*A + MB.*C;
        MBT = MA.*B + MB.*D;
        MCT = MC.*A + MD.*C;
        MDT = MC.*B + MD.*D;
        MA = MAT; MB = MBT; MC = MCT; MD = MDT;
      end
      nHole = nHole + 1;
    end
    
    % Cascade cylindrical or conical sections
    if L(m) < eps, continue; end % skip if at a diameter discontinuity
    if ra(m) == ra(m+1)
      [Gamma, ZcLoss] = lossesCylinder(k, ra(n), Zc(n), c, rho, gamma, lv, Pr, lossy, alphacm);
      [A, B, C, D] = tmmCylinder( Gamma, L(n), ra(n), ZcLoss);
    else
      [Gamma, ZcLoss] = lossesCone(k, ra(n),ra(n+1), L(n), Zc(n), c, rho, gamma, lv, Pr, lossy, alphacm);
      [A, B, C, D] = tmmCone(k, Gamma, L(n), ra(n), ra(n+1), ZcLoss);
    end
    MAT = MA.*A + MB.*C;
    MBT = MA.*B + MB.*D;
    MCT = MC.*A + MD.*C;
    MDT = MC.*B + MD.*D;
    MA = MAT; MB = MBT; MC = MCT; MD = MDT;
  end

  % Finally, include 1/2 of series length correction of next open hole
  % (if not the end hole).
  if n < nOth
    nHole = oidx(n+1);
    [~, B, ~, ~] = tmmTonehole( k, rb(nHole)/ras(nHole), rb(nHole), ...
      t(nHole), Zch(nHole), states(nHole), wallcst, alphacm, '', ...
      chimney(nHole), padr(nHole), padt(nHole), holew(nHole) );
    B = B/2; C = 0; A = 1; D = 1;
    MAT = MA.*A + MB.*C;
    MBT = MA.*B + MB.*D;
    MCT = MC.*A + MD.*C;
    MDT = MC.*B + MD.*D;
    MA = MAT; MB = MBT; MC = MCT; MD = MDT;
  end

  % Fill Y matrix from cascaded terms.
  if n > 1
    Y(n,n-1,:) = Ymunm1;
  end
  Y(n,n,:) = Ypnm1 + MD./MB;
  if n < nOpen
    Y(n,n+1,:) = -1./MB;
    Ymunm1 = -1./MB;
    Ypnm1 = MA./MB;
  end

end

% Now handle end condition
if endType
  switch endType
    case 1
      ZB(nOpen,nOpen,:) = Zc(end)*radiation( k, ra(end), 'UnflangedDalmont'); % self-radiation
    case 2
      ZB(nOpen,nOpen,:) = Zc(end)*radiation( k, ra(end), 'FlangedNorris'); % self-radiation
    case 3
      ZB(nOpen,nOpen,:) = zeros(size(k));
  end
  if doInteractions
    for m = 1:nOth
      dnm = abs( x(end) - ds(m) );
      ZB(nOpen, m, :) = 1j*rho*c*k.*exp(-1j*k*dnm)/(4*pi*dnm);
    end
  end
  Y(nOpen,nOpen-1,:) = Ymunm1;
  Y(nOpen,nOpen,:) = Ypnm1;
else % end is closed
  Y(nOth,nOth,:) = Ypnm1 + MC./MA;
end

% Compute Eqs. 14 & 10 for each frequency
I = eye(nOpen);
U = zeros(nOpen, length(k));
P = zeros(nOpen, length(k));
Us = zeros(nOpen, 1);
Us(1) = 1;
for n = 1:length(k)
  U(:, n) = ( I + Y(:, :, n)*ZB(:, :, n) )^(-1) * Us;
  P(:, n) = ZB(:, :, n) * U(:, n);
end
  
% Compute "load" impedance of section rightward from first open hole
%Zl = (P(1, :) ./ U(1, :));
Zl = P(1, :);  % since Us(1) = 1

% If there is at least one open tonehole, account for 1/2 of its series
% impedance on the upstream side.
if ~isempty( oidx )
  nHole = oidx(1);
  [~, B, ~, D,] = tmmTonehole( k, rb(nHole)/ras(nHole), rb(nHole), ...
    t(nHole), Zch(nHole), states(nHole), wallcst, alphacm, '', ...
    chimney(nHole), padr(nHole), padt(nHole), holew(nHole) );
  Zl = (Zl + B) ./ D;
  nHole = nHole - 1; % decrement hole counter
else % no open holes
  nHole = sum(isHole);
end

% Continue back through the remainder of the closed system using the TMM
% (from the first closed hole upstream of most upstream open hole).
for n = xidx(1)-1:-1:1
  if L(n) > eps
    if ( ra(n) == ra(n+1) )
      [A, B, C, D] = tmmCylinder( k, L(n), ra(n), Zc(n), wallcst, alphacm );
    else
      [A, B, C, D] = tmmCone( k, L(n), ra(n), ra(n+1), Zc(n), wallcst, alphacm );
    end
    Zl = (A.*Zl + B) ./ (C.*Zl + D);
  end
  
  if isHole(n)
    [A, B, C, D] = tmmTonehole( k, rb(nHole)/ra(n), rb(nHole), t(nHole), ...
      Zch(nHole), states(nHole), wallcst, alphacm, '', chimney(nHole), ...
      padr(nHole), padt(nHole), holew(nHole) );
    nHole = nHole - 1;
    Zl = (A.*Zl + B) ./ (C.*Zl + D);
  end
end

Zin = Zl / Zc(1);
