clear all; clc;
load('modele_transmitancje.mat');

figure('Name', 'Weryfikacja modeli', 'Color', 'w');
hold on; grid on;
title('Weryfikacja identyfikacji - porównanie odpowiedzi obiektu i modeli');
xlabel('Czas t [s]');
ylabel('Temperatura wewnętrzna T_{wew} [°C]');

% Rysowanie czystych przebiegów odpowiedzi
h1 = plot(t, aTwew, 'b-', 'LineWidth', 2.5);
h2 = plot(t, Twew_model_2p, 'm-.', 'LineWidth', 2.0);
h3 = plot(t, Twew_model_st, 'g--', 'LineWidth', 2.0);
h4 = plot(t, Twew_model_mom, 'r:', 'LineWidth', 2.0);

ylim([y0 - 0.05, y_ust + 0.1]);
legend([h1, h2, h3, h4], ...
       {'Obiekt rzeczywisty (T_{wew})', ...
        sprintf('Model dwupunktowy (K=%.1f, T=%.1f, T_0=%.1f)', K_2p, T_2p, T0_2p), ...
        sprintf('Model ze stycznej (K=%.1f, T=%.1f, T_0=%.1f)', K_st, T_st, T0_st), ...
        'Model z metody momentów (II rzędu)'}, ...
       'FontSize', 10, 'Location', 'southeast');