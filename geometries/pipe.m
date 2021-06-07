function [boreData, holeData] = pipe( ~ )
% Air column and tonehole data for a straight cylinderical pipe with no
% side holes.

% Air column radii measurements
boreData = [1e-3*[0     599];   % positions (from input end)
            1e-3*[15.5  15.5]/2];   % radii at corresponding points

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

holeData = zeros(6, 0);
