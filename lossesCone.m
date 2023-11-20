function [G, Zc] = lossesCone(k, R1, R2, L, Zc, c, rho, gamma, lv, Pr, lossy, alphacm)
% [G, ZC] = LOSSESCONE( K, R1, R2, L, ZC, C, RHO, GAMMA, LV, PR, LOSSY,
% ALPHACM) computes the wave propagation variable GAMMA and the
% characteristic impedance ZC for given values of K = omega/c in a
% conical section with input radius R1 and output radius R2 and length L. 
% C, RHO, GAMMA, LV, and PR are physical constants corresponding to the
% speed of sound, density of air, ratio of specific heats, viscous boundary
% layer thickness, and Prandtl number. ALPHACM is an optional 1D vector 
% additional loss factor (the same size as K) of attenuation values
% corresponding to molecular and classical losses in air. 
% see PHYSICALSETTINGS for more details. If parameter LOSSY = 0, losses 
% will be ignored. Otherwise various approximations, with increasing 
% accuracy, can be used instead:
%   1 - "Standard" TMM losses from [1], [2], [3]. Zeroth/First order
%   approximations. 
%   2 - Second/Fourth order approximations (depending on bore radius) 
%   from [4], [5] 
%   3 - Numerical computation of Bessel functions of thie first kind, from
%   [4]
% For LOSSY = 1 and LOSSY = 2. The type of loss (wide or narrow pipe) will
% change based on the maximum Stokes number, rv = R*sqrt(k/lv). If max(rv)
% < 1 then the narrow pipe approximation will be used instead
%
% NOTE: Currently, this function uses approximations for cylindrical losses 
% with an assumed equivalent radius [1], [5]. New research has been
% published on more accurate conical section losses. This should be
% implemented instead. 
% 
% by Champ Darabundit, McGill University, 2023.
%
% References:
%
% 1. A. Lefebvre, G.P. Scavone, and J. Kergomard (2013), "External tonehole
%       interactions in woodwind instruments." Acta Acustica united with
%       Acustica, vol. 99, pp. 975-985
%
% 2. J. Abel, T. Smyth, and J.O. Smith (2003), "A simple, accurate wall
%       loss filter for acoustic tubes." in proc. of the 6th Int. Conf. 
%       Digital Audio Effects (DAFx-03), London, UK. 
%
% 3. A.H. Benade (1968), "On the propagation of sound waves in a
%       cylindrical conduit." Journal of the Acoustical Society of America,
%       vol. 44, pp. 616-623. 
%
% 4. C. Zwikker and C.W. Kosten (1949), "Sound absorbing materials."
%       Elsevier Publishing Company. 
%
% 5. A. Chaigne and J. Kergomard (2016), "Acoustics of music instruments."
%       Springer-Verlag. 

if lossy > 3 || mod(lossy,1) ~= 0
    error('lossy must be an integer between 0 and 3')
end

% Compute equivalent radius
x1 = R1 * L / (R2 - R1);
Req = L * R1 / x1 / log(1 + L/x1);

% Case 0: No losses
if lossy == 0
    G = 1j*k;
    Zc = R1 * Zc / R2;  % Scale Zc

% Case 1: TMM approximation from [1], [2], [3]
elseif lossy == 1
    rv = Req*sqrt(abs(k/lv));
    Zc = R1 * Zc / R2;  % Scale Zc
    if max(rv) < 1
        % G = k.*2.*sqrt(gamma)./rv + 1j.*k./rv.*2.*sqrt(gamma);
        G = k.*2.*sqrt(gamma)./rv.*(1 + 1j);
    else
        cst = sqrt(lv/2).*(1 + (gamma - 1)/sqrt(Pr));
        G = 1j*k + (1 + 1j) * cst .* sqrt(k)/Req;
    end

% Case 2: 2nd/4th order expansion approximations from [4], [5]
elseif lossy == 2
    lt = lv/Pr;
    rv = Req*sqrt(abs(k/lv));
    rt = Req*sqrt(abs(k/lt));
    S = pi*Req*Req;
    w = k.*c;
    % Narrow pipe, this condition could be changed...
    if max(rv) < 1  
        Zv = (rho*c/S).*(8*lv)./Req^2.*(1 + (1/6).*1j.*rv.^2);
        Yt = 1j.*w.*gamma.*(S./(rho*c^2)).*(1 - (gamma - 1)./(8*gamma).*1j.*rt.^2);
    % Wide pipe
    else            
        Zv = (1j*w*rho./S).*(1 + (2/sqrt(2))*(1 - 1j)./rv - 3j./rv.^2);
        Yt = (1j*w*S)./(rho*c^2).*(1 + (gamma - 1).*((2/sqrt(2))*(1 - 1j)./rt + 1j./rt.^2));
    end
    G = sqrt(Zv.*Yt);
    Zc = sqrt(Zv./Yt);

% Case 3: Numerical computation of Bessel function from [4], [5]
elseif lossy == 3
    lt = lv/Pr;
    kv = ((1 - 1j)/sqrt(2))*sqrt(k/lv);
    kt = ((1 - 1j)/sqrt(2))*sqrt(k/lt);
    S = pi*Req*Req;
    w = k.*c;
    % Compute Bessel functions with scaling in case of overflow
    J1v = besselj(1, kv.*Req, 1); J1v = J1v./ exp(-abs(imag(kv.*Req)));
    J0v = besselj(0, kv.*Req, 1); J0v = J0v./ exp(-abs(imag(kv.*Req)));
    J1t = besselj(1, kt.*Req, 1); J1t = J1t./ exp(-abs(imag(kt.*Req)));
    J0t = besselj(0, kt.*Req, 1); J0t = J0t./ exp(-abs(imag(kt.*Req)));
    % Compute immittances
    Zv = (1j.*w.*rho/S).*(1 - 2./(kv.*Req).*(J1v./J0v)).^(-1);
    Yt = (1j.*w.*S)/(rho.*c^2).*(1 + (gamma - 1).*(2./(kt.*Req)).*(J1t./J0t));
    G = sqrt(Zv.*Yt);
    Zc = sqrt(Zv./Yt);
end

% Include alphacm if used
if exist( 'alphacm', 'var')
  if size(k) ~= size(alphacm)
    error( 'Incompatible sizes of k and alphacm vectors.' );
  end
  G = G + alphacm;
end

end