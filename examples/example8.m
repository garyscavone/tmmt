% A matlab script used to compare the lossy wavenumber Gamma and
% characteristic impedance Zc for different bore radii. Useful to see range
% of validity for different approximations
%
% by Champ Darabundit, McGill University, 2023

clear all; close all; clc; warning('off')

% Evaluation frequencies
N = 2^11;                   % Number of evaluation points
fs = 48e3;                  % Sampling frequency, evaluate up to Nyquist
f = [0:N/2]*(fs/N);

T = 20;                     % temperature (C)

% Include path needed to script
addpath('../')

% Get physical constants based on temperature
[c, rho, gamma, lv, Pr, alphacm] = physicalSettings(T, f);
k = (2*pi*f)./c;            % Compute the wavenumber

% Radii to compare and initial characteristic impedance
ra = [5e-3, 5e-4, 5e-5, 5e-6];
Zc = rho.*c./(pi.*ra.^2);
% Linestyle for plot
ls = {'--', ':', '-'};
for n = 1:4
    % Keep track of viscous boundary-layer ratio to see range of validity
    % for different approximations
    rv = ra(n)*sqrt(abs(k/lv));
    figure()
    for lossy = 3:-1:1
        % Do loss calculations
        [G, ZcLoss] = lossesCylinder(k, ra(n), Zc(n), c, rho, gamma, lv, Pr, lossy);

        % If Zc is scalar, multiply by ones for the plot
        if numel(ZcLoss) == 1
            ZcLoss = ones(size(f)).*ZcLoss;
        end

        % Plot and compare different loss approximations
        subplot(221)
        plot(f, real(G), LineStyle=ls{lossy}, LineWidth=2)
        hold on
        legend('3: Bessel.', '2: CK Approx.', '1: TMM Approx.', Location='best')
        grid on
        xlim([0, fs/2])
        xlabel('Frequency (Hz)')
        ylabel('\Re\{\Gamma\}')
        title('Lossy Wavenumer \Gamma')

        subplot(223)
        plot(f, imag(G), LineStyle=ls{lossy}, LineWidth=2)
        hold on
        legend('3: Bessel.', '2: CK Approx.', '1: TMM Approx.', Location='best')
        grid on
        xlim([0, fs/2])
        xlabel('Frequency (Hz)')
        ylabel('\Im\{\Gamma\}')
        title('Lossy Wavenumer \Gamma')

        subplot(222)
        plot(f, real(ZcLoss), LineStyle=ls{lossy}, LineWidth=2)
        hold on
        legend('3: Bessel.', '2: CK Approx.', '1: TMM Approx.', Location='best')
        grid on
        xlim([0, fs/2])
        xlabel('Frequency (Hz)')
        ylabel('\Re\{Z_c\}')
        title('Char. Impedance Z_c')

        subplot(224)
        plot(f, imag(ZcLoss), LineStyle=ls{lossy}, LineWidth=2)
        hold on
        legend('3: Bessel.', '2: CK Approx.', '1: TMM Approx.', Location='best')
        grid on
        xlim([0, fs/2])
        xlabel('Frequency (Hz)')
        ylabel('\Im\{Z_c\}')
        title('Char. Impedance Z_c')
    end
    % Title for each radii comparison
    sgtitle("Loss Comparison: " + sprintf("R = 5e{-%i} m", n+2) + ", r_v \in " + sprintf("[%.2f, %.2f]", min(rv(rv > 0)), max(rv)))
    saveas(gcf, sprintf("5e-%i.png", n+2))
end
