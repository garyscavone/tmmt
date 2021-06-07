function [boreData, holeData] = keefeFlute( fingering )
% Air column, tonehole, and fingering data for a flute given
% by Doug Keefe in "Woodwind air column models," JASA, 1990.  The
% instrument has 6 toneholes which are cut directly into the bore
% without pads or keys.

if ( nargin < 1 || fingering < 1 || fingering > 7 )

  fingerInfo = [
    '                                  ';
    'fingering 1: [ 0 0 0 0 0 0 ]; % D ';
    'fingering 2: [ 0 0 0 0 0 1 ]; % E ';
    'fingering 3: [ 0 0 0 0 1 1 ]; % F ';
    'fingering 4: [ 0 0 0 1 1 1 ]; % G ';
    'fingering 5: [ 0 0 1 1 1 1 ]; % A ';
    'fingering 6: [ 0 1 1 1 1 1 ]; % B ';
    'fingering 7: [ 1 1 1 1 1 1 ]; % C ';
    '(closed state = 0, open state = 1)';
    '                                  ';
    ];

  disp( fingerInfo );
  fingering = 0;
  while ( fingering < 1 || fingering > 7 )
    choice = inputdlg('Select a fingering: ', 'Fingering', 1, {num2str(1)} );
    if isempty( choice )
      boreData = []; holeData = [];
      return;
    end
    fingering = str2double( choice{1} );
  end
  
end

% Air column radius measurements
boreData = [1e-3*[0     575.2];     % positions (from input end)
            1e-3*[18.9  18.9]/2];   % radii at corresponding points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tonehole specifications
%
% holeData: Row 1 - positions (from input end)
%           Row 2 - tonehole radii
%           Row 3 - tonehole heights
%           Row 4 - tonehole protrusion length out of bore
%           Row 5 - tonehole states (1 = open, 0 = closed)
%           Row 6 - tonehole pad state (1 = pad exists, 0 = no pad
%                ** subsequent rows optional **
%           Row 7 - tonehole pad radii
%           Row 8 - tonehole pad heights (above hole)
%           Row 9 - tonehole wall thickness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

holeData = [1e-3*[286.4 323.4 359  412  436.4 475.7];
            1e-3*[9.53  9.53  7.94 7.94 9.53  6.35 ]/2;
            1e-3*[3.4   3.4   3.4  3.4  3.4   3.4  ];
                 [0     0     0    0    0     0    ];
                 [0     0     0    0    0     0    ];
                 [0     0     0    0    0     0    ]];

% Fingerings: 1 = tonehole open, 0 = tonehole closed
fingerings(1,:) = [ 0 0 0 0 0 0 ]; % D
fingerings(2,:) = [ 0 0 0 0 0 1 ]; % E
fingerings(3,:) = [ 0 0 0 0 1 1 ]; % F
fingerings(4,:) = [ 0 0 0 1 1 1 ]; % G
fingerings(5,:) = [ 0 0 1 1 1 1 ]; % A
fingerings(6,:) = [ 0 1 1 1 1 1 ]; % B
fingerings(7,:) = [ 1 1 1 1 1 1 ]; % C

holeData(5,:) = fingerings( fingering, : );

