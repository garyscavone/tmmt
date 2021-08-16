function rzplot( f, Zin, plotTypes, doGrid, doHold, xlimits, pcolor, c, doSmooth )
% RZPLOT: Plot impedance and reflectance (and time-domain correlates).
%
% RZPLOT( f, Zin, plotTypes, doGrid, doHold, xlimit, pcolor, c, doSmooth ) where:
%   - f: 1D vector of frequencies;
%   - Zin: 1D vector of impedance data;
%   - plotTypes: 1 or more values of the following plot types:
%     - 1:  Impedance Magnitude in dB
%     - 2:  Impedance Magnitude (logarithmic y-scale)
%     - 3:  Impedance Phase
%     - 4:  Impedance Real Part
%     - 5:  Impedance Imaginary Part
%     - 6:  Reflectance Magnitude
%     - 7:  Reflectance Phase
%     - 8:  Reflectance Real Part
%     - 9:  Reflectance Imaginary Part
%     - 10: Reflectance Equivalent Length
%     - 11: Impulse Response
%     - 12: Reflection Function
%   - doGrid: optional boolian flag to turn on the grid;
%   - doHold: optional boolian flag to specify "hold on" before plotting;
%   - xlimits: optional x-axis plot limits;
%   - pcolor: optional color for plot data;
%   - c: optional speed of sound value needed for equivalent length plot;
%   - doSmooth: optional boolian flag to specify smoothing of time data.
%
% by Gary P. Scavone, McGill University, 2013-2021.

if size( Zin) ~= size( f )
  error( 'Zin and f must be the same size.' );
end

if ~isvector( Zin )
  error( 'Zin and f must be a 1D vectors.' );
end

if ( sum(plotTypes > 4) )
  R = (Zin - 1) ./ (Zin + 1);
end

if ~exist( 'doHold', 'var' )
  doHold = false;
end

if ~exist( 'pcolor', 'var' )
  pcolor = 'b';
end

if ~exist( 'doGrid', 'var' )
  doGrid = false;
end

if ~exist( 'doSmooth', 'var' )
  doSmooth = false;
end

nPlots = length( plotTypes );
for n = 1:nPlots
  subplot(nPlots, 1, n)
  if doHold, hold on; end
  if plotTypes(n) == 1
    p = plot( f, 20*log10(abs(Zin)), pcolor );
    ylabel('20*log10(|Impedance|)')
  elseif plotTypes(n) == 2
    p = semilogy( f, abs(Zin), pcolor);
    ylim([0.01 100])
    ylabel('|Impedance|')
  elseif plotTypes(n) == 3
    p = plot( f, angle(Zin), pcolor );
    ylabel('Impedance Phase')
  elseif plotTypes(n) == 4
    p = plot( f, real(Zin), pcolor );
    ylabel('Real(Impedance)')
  elseif plotTypes(n) == 5
    p = plot( f, imag(Zin), pcolor );
    ylabel('Imag(Impedance)')
  elseif plotTypes(n) == 6
    p = plot( f, abs(R), pcolor );
    ylabel('Reflectance')
    ylabel('|R|')
  elseif plotTypes(n) == 7
    p = plot( f, angle(R), pcolor );
    ylabel('Reflectance Phase')
  elseif plotTypes(n) == 8
    p = plot( f, real(R), pcolor );
    ylabel('Real(Reflectance)')
  elseif plotTypes(n) == 9
    p = plot( f, imag(R), pcolor );
    ylabel('Imag(Reflectance)')
  elseif plotTypes(n) == 10
    if ~exist( 'c', 'var' )
      c = 340;
    end
    el = (pi - unwrap(angle(R)) + 1j*log(abs(R))) ./ (4*pi*f/c);
    p = plot( f, 1000*(real(el)), pcolor );
    ylim([150 200])
    ylabel('Equivalent Length (mm)')
  end
  if ( plotTypes(n) < 11 )
    hx = xlabel('Frequency (Hz)');
  else
    if plotTypes(n) == 11
      H = transpose(Zin);
      H(1) = 0;           % impedance at f=0 should be zero (roughly)
      ytext = 'Impulse Response';
    elseif plotTypes(n) == 12
      H = transpose(R);
      H(1) = 1;           % reflectance at f=0 should be one (roughly)
      ytext = 'Reflection Function';
    end
    H(end) = real(H(end));
    if isrow( H )
      H = [ H, conj( flip( H(2:end-1) ) ) ];  % make conjugate symmetric
    else
      H = [ H; conj( flip( H(2:end-1) ) ) ];  % make conjugate symmetric
    end
    h = real( ifft( H ) );
    t = (0:length(h)-1) * 500 / max(f);      % time in milliseconds %
    if doSmooth
      h = filtfilt([0.5 0.5], 1, h);
    end
    p = plot( t, h, pcolor );
    hx = xlabel('Time (ms)');
    ylabel( ytext );
    hold off;
  end
  if doGrid
    grid on;
  end
  if exist( 'xlimits', 'var' ) && ~isempty( xlimits )
    xlim( xlimits );
  end
end

%set(p, 'LineWidth', 1.5)
%set( [gca, hx], 'fontsize', 18, 'fontname', 'Times' );
%set(p,'MarkerSize',1.5)