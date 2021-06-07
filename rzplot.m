function rzplot( f, Zin, plotTypes, doGrid, doHold, xlimits, pcolor, c, doSmooth )
% RZPLOT: Plot impedance and reflectance (and time-domain correlates).
%
% RZPLOT( f, Zin, plotTypes, doGrid, doHold, xlimit, pcolor, c, doSmooth ) where:
%   - f: 1D vector of frequencies;
%   - Zin: 1D vector of impedance data;
%   - plotTypes: 1 or more values of the following plot types:
%     - 1:  Impedance Magnitude in dB
%     - 2:  Impedance Phase
%     - 3:  Impedance Real Part
%     - 4:  Impedance Imaginary Part
%     - 5:  Reflectance Magnitude
%     - 6:  Reflectance Phase
%     - 7:  Reflectance Real Part
%     - 8:  Reflectance Imaginary Part
%     - 9:  Reflectance Equivalent Length
%     - 10: Impedance Magnitude (logarithmic y-scale)
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

[M, N] = size( Zin );
if ( M == 1 )
  Zin = Zin.';
elseif ( M > 1 && N > 1 )
  error( 'Zin should be a 1D vector.' );
end

[M, N] = size( f );
if ( M == 1 )
  f = f.';
elseif ( M > 1 && N > 1 )
  error( 'f should be a 1D vector.' );
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
    p = plot( f, angle(Zin), pcolor );
    ylabel('Input Impedance Phase')
  elseif plotTypes(n) == 3
    p = plot( f, real(Zin), pcolor );
    ylabel('Real(Input Impedance)')
  elseif plotTypes(n) == 4
    p = plot( f, imag(Zin), pcolor );
    ylabel('Imag(Input Impedance)')
  elseif plotTypes(n) == 5
    p = plot( f, abs(R), pcolor );
    ylabel('Reflectance')
    ylabel('|R|')
  elseif plotTypes(n) == 6
    p = plot( f, angle(R), pcolor );
    ylabel('Reflectance Phase')
  elseif plotTypes(n) == 7
    p = plot( f, real(R), pcolor );
    ylabel('Real(Reflectance)')
  elseif plotTypes(n) == 8
    p = plot( f, imag(R), pcolor );
    ylabel('Imag(Reflectance)')
  elseif plotTypes(n) == 9
    if ~exist( 'c', 'var' )
      c = 340;
    end
    el = (pi - unwrap(angle(R)) + 1j*log(abs(R))) ./ (4*pi*f/c);
    p = plot( f, 1000*(real(el)), pcolor );
    ylim([150 200])
    ylabel('Equivalent Length (mm)')
  elseif plotTypes(n) == 10
    p = semilogy( f, abs(Zin), pcolor);
    ylim([0.01 100])
    ylabel('|Input Impedance| / Z0')
    ylabel('| Z / Z_0 |')
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
    H = [ H, conj( fliplr( H(2:end-1) ) ) ];  % make conjugate symmetric
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

set(p,'LineWidth',1.5)
%set( [gca, hx], 'fontsize', 18, 'fontname', 'Times' );
%set(p,'MarkerSize',1.5)