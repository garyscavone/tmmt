function drawBore( instrument, fingering, boreData, holeData )
% DRAWBORE:  Make a 3D plot of the instrument specified by the input
%            arguments.
%
% DRAWBORE( INSTRUMENT, FINGERING, BOREDATA, HOLEDATA ) draws a 3D plot of
% the instrument geometry specified either in the file INSTRUMENT or in the
% BOREDATA and HOLEDATA vectors (when INSTRUMENT is an empty array). The
% FINGERING is optional (default = 1). INSTRUMENT must be a string.
%
% by Gary P. Scavone, McGill University, 2013-2020.

if ( nargin < 1 )
  [instrument, ~] = uigetfile('*.m', 'Choose an instrument file:');
  if instrument == 0
    return;
  end
  instrument = instrument(1:strfind(instrument, '.m')-1);
end

if ( nargin < 2 )
  fingering = 1;
end

if ~isempty( instrument )
  [boreData, holeData] = feval( instrument, fingering );
else
  if nargin < 3
    error('boreData vector needed when no instrument specified.');
  end
  if nargin < 4
    holeData = zeros(6, 0);
  end
end

xMax = max( boreData(1, :) );
rMax = max( boreData(2,:) );

clf
axis([0 xMax -2*rMax 2*rMax -2*rMax 2*rMax]);
view(-10, 55);
colormap('bone');
hold on;
grid on;
xlabel('x-axis of instrument');
zlabel('z-axis of instrument');

% Draw main air column
nSections = length(boreData(1, :));
L = diff( boreData(1, :) );
[zb, yb, xb] = cylinder( boreData(2, 1:2) );
zb = zb(1, :);
yb = yb(1, :);
xb = xb(1, :);
for n = 2:nSections
  [z, y, x] = cylinder( boreData(2, n-1:n) );
  zb = [zb; z(2, :)];
  yb = [yb; y(2, :)];
  xb = [xb; xb(n-1, 1)+(x(2, :)*L(n-1))];
end
surf(xb, yb, zb);

% Draw toneholes
idx = find(diff(boreData(1,:)) == 0);
boreData(1, idx+1) = boreData(1, idx+1) + eps; % avoid double values
lastx = 0;
counter = 0;
for n = 1:length(holeData(1,:))
  [xh, yh, zh] = cylinder(holeData(2,n));
  ra = interp1(boreData(1,:), boreData(2,:), holeData(1,n), 'linear');
  t = holeData(3, n);
  if holeData(4,n) == 0 % closed - black
    fillcolor = [0 0 0];
  else
    fillcolor = [1 1 1]; % open - white
  end
  % Avoid overlapping holes
  if lastx > (holeData(1,n) - holeData(2,n))
    counter = counter + 1;
    dscale = -1;
    if ( mod( counter, 2 ) == 0 ) dscale = 1; end
    mesh(xh+holeData(1,n), dscale*(ra+zh*t), yh);
    fill3(xh(2,:)+holeData(1,n), dscale*(ra+zh(2,:)*t), yh(2,:), fillcolor);
  else
    counter = 0;
    mesh(xh+holeData(1,n), yh, ra+zh*t);
    fill3(xh(2,:)+holeData(1,n), yh(2,:), ra+zh(2,:)*t, fillcolor);
    lastx = holeData(1,n) + holeData(2,n);
  end
end

hold off