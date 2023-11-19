function [G, Zc] = lossesCylinder(k, R, Zc, c, rho, gamma, lv, Pr, lossy, alphacm)

if lossy == 0
    G = 1j*k;
elseif lossy == 1
    cst = sqrt(lv/2).*(1 + (gamma - 1)/sqrt(Pr));
    G = 1j*k + (1 + 1j) * cst .* sqrt(k)/R;
elseif lossy == 2
    lt = lv/Pr;
    rv = R*sqrt(abs(k/lv));
    rt = R*sqrt(abs(k/lt));
    S = pi*R*R;
    w = k.*c;
    % Narrow pipe, this condition could be changed...
    if max(rv) < 1  
        Zv = (rho*c/S).*(8*lv)./R^2.*(1 + (1/6).*1j.*rv.^2);
        Yt = 1j.*w.*(S./(rho*c^2)).*(1 - (gamma - 1)./(8*gamma).*1j.*rt.^2);
    % Wide pipe
    else            
        Zv = (1j*w*rho./S).*(1 + (2/sqrt(2))*(1 - 1j)./rv - 3j./rv.^2);
        Yt = (1j*w*S)./(rho*c^2).*(1 + (gamma - 1).*((2/sqrt(2))*(1 - 1j)./rt + 1j./rt.^2));
    end
    G = sqrt(Zv.*Yt);
    Zc = sqrt(Zv./Yt);
elseif lossy == 3
    lt = lv/Pr;
    kv = ((1 - 1j)/sqrt(2))*sqrt(k/lv);
    kt = ((1 - 1j)/sqrt(2))*sqrt(k/lt);
    S = pi*R*R;
    w = k.*c;
    % Compute Bessel functions with scaling in case of overflow
    J1v = besselj(1, kv.*R, 1); J1v = J1v./ exp(-abs(imag(kv.*R)));
    J0v = besselj(0, kv.*R, 1); J0v = J0v./ exp(-abs(imag(kv.*R)));
    J1t = besselj(1, kt.*R, 1); J1t = J1t./ exp(-abs(imag(kt.*R)));
    J0t = besselj(0, kt.*R, 1); J0t = J0t./ exp(-abs(imag(kt.*R)));
    % Compute immittances
    Zv = (1j.*w.*rho/S).*(1 - 2./(kv.*R).*(J1v./J0v)).^(-1);
    Yt = (1j.*w.*S)/(rho.*c^2).*(1 + (gamma - 1).*(2./(kt.*R)).*(J1t./J0t));
    G = sqrt(Zv.*Yt);
    Zc = sqrt(Zv./Yt);
end

if exist( 'alphacm', 'var')
  if size(k) ~= size(alphacm)
    error( 'Incompatible sizes of k and alphacm vectors.' );
  end
  G = G + alphacm;
end

end