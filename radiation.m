function Zr = radiation( k, a, type )
% ZR = RADIATION( K, A, TYPE ) computes the radiation impedance (normalized
%      by Zc) for the given values of k = omega / c and a = output radius
%      in meters. TYPE is an optional parameter specifying a particular
%      condition or formula. The default is 'UnflangedDalmont,' which is an
%      approximation provided in [1].
%
% by Gary P. Scavone, McGill University, 2013-2021.
% Based in part on functions from WIAT by Antoine Lefebvre.
%
% References:
%
% 1. J. Dalmont and C.J. Nederveen, "Radiation impedance of tubes with
%      different flanges: numerical and experimental investigations,"
%      Journal of Sound and Vibration,  Vol. 244, pp. 505-534, 2001.
%
% 2. R. Caussé, J. Kergomard, and X. Lurton, "Input impedance of brass
%      musical instruments - Comparaison between experimental and numerical
%      models," J. Acoust. Soc. Am.,  Vol. 75, pp. 241-254, 1984.
%
% 3. H. Levine and J. Schwinger, "On the radiation of sound from an
%      unflanged circular pipe," Phys. Rev., 73(4), pp. 383-406, 1948.
%
% 4. A.N. Norris and I.C. Sheng. "Acoustic radiation from a circular pipe
%      with an infinite flange." Journal of Sound and Vibration, Vol. 135,
%      pp. 85-93, 1989.
%
% Generally, Zr = Zc*(1+R)/(1-R), where R is the reflection coefficient for
% the open end of a pipe and Zc is the characteristic impedance.  R can
% also be expressed as:
%     R = -|R0|exp(-2j*k*a*delta),
% where delta is the frequency dependant length correction.  See
% reference [1], pg. 509.  A low-frequency approximation for an unflanged
% Zr is:
%     Zr = 0.25*ka^2 + 0.61j*ka

if ( nargin < 3 )
  type = 'UnflangedDalmont';
end

ka = k*a;
ka2 = ka.^2;

if strcmp( type, 'unflanged' )
  % The Levine & Schwinger results are calculated by numerical
  % integrations. Eq. (VI.5) of Levine and Schwinger is used to calculate
  % the reflection coefficient magnitude.
  intvi5 = @(x, z) atan(besselk(1, x)./(pi*besseli(1, x))) ...
    .*(1 - z./sqrt(z^2 + x.^2)).*(1./x);

  if ka(1) == 0, ka(1) = eps; end % avoid division by 0
  sum1 = zeros(size(ka));
  for n=1:length(ka)
    sum1(n) = integral(@(x)intvi5(x, ka(n)), 0, 20); % upper limit determined empirically
  end
  r = (pi*ka).^(0.5).*exp(-ka+(sum1./pi));
  
  % First and second integrals in Eq. (VI.4) of Levine and Schwinger, used
  % to calculate the length correction.
  intvi4a = @(x, z) log(pi*abs(besselj(1,x)).*sqrt(besselj(1,x).^2 + ...
    bessely(1,x).^2))./(x.*sqrt(z^2 - x.^2));
  intvi4b = @(x, z) log(1./(2*besseli(1,x).*besselk(1,x))) ...
    ./(x.*sqrt(x.^2 + z^2));

  warning('off'); % turn off warning about integration interval if ka(1) close to zero.
  sum2 = zeros(size(ka));
  for n=1:length(ka)
    sum2(n) = integral(@(x)intvi4a(x, ka(n)), 0, ka(n));
  end
  sum3 = zeros(size(ka));
  for n=1:length(ka)
    sum3(n) = integral(@(x)intvi4b(x, ka(n)), 0, 600); % upper limit determined empirically
  end
  warning('on')
  loa = (sum2+sum3)./pi;
  R = -r.*exp(-2*1i*ka.*loa);
  Zr = (1 + R) ./ (1 - R);

elseif strcmp( type, 'UnflangedDalmont' )
  
  % Unflanged pipe radiation impedance approximation (ka < 1.5) from
  % reference [1].
  R0 = (1 + 0.2*ka - 0.084*ka2) ./ (1 + 0.2*ka + (0.5-0.084)*ka2); % ka<3.5
  delta = 0.6133*((1 + 0.044*ka2)./(1 + 0.19*ka2) - 0.02*sin(2*ka).^2);
  R = -abs(R0).*exp(-2*1i*ka.*delta);
  Zr = (1 + R) ./ (1 - R);

elseif strcmp( type, 'UnflangedCausse' )
  
  % Unflanged pipe radiation impedance approximation (ka < 1.5) from
  % reference [2], formula taken from reference [1].
  Zr = 1j*0.6113*ka - 1j*ka.^3 .* (0.036-0.034*log(ka) + 0.0187*ka2) + ...
    0.25*ka2 + ka.^4.*(0.0127+0.082*log(ka) - 0.023*ka2);
  
elseif strcmp( type, 'FlangedNorris' )

  % Infinite flange fit formula by Norris and Sheng for ka < 3.5, taken
  % from reference [1].
  %
  % A.N. Norris and I.C. Sheng. "Acoustic radiation from a circular pipe
  % with an infinite flange." Journal of Sound and Vibration, Vol. 135,
  % pp. 85-93, 1989.
  R0 = (1 + 0.323*ka - 0.077*ka2) ./ (1 + 0.323*ka + (1-0.077)*ka2);
  delta = 0.8216*(1 + (0.77*ka).^2./(1 + 0.77*ka)).^(-1);
  R = -abs(R0).*exp(-2*1i*ka.*delta);
  Zr = (1 + R) ./ (1 - R);

else
  error( 'Unknown type argument.' );
end