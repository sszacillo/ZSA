% =========================================================================
% SKRYPT GENERUJĄCY I EKSPORTUJĄCY WYKRESY DO LATEXA (300 DPI, BEZ RAMEK)
% =========================================================================
clear all; close all; clc;

load('modele_transmitancje.mat');

% --- 1. Wykres: metoda_dwupunktowa.png ---
f1 = figure('Color', 'w', 'Position', [100, 100, 950, 480]); 
hold on; grid on; box on;
title('Identyfikacja metodą dwupunktową - konstrukcja graficzna', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11); 
ylabel('Temperatura wewnętrzna T_{wew} [°C]', 'FontSize', 11);

h_ob = plot(t, aTwew, 'b-', 'LineWidth', 2.2);
yline(y0, 'k--', 'Poziom początkowy T_{wew0}', 'LabelHorizontalAlignment', 'left');
yline(y_ust, 'k--', 'Wartość ustalona T_{wew,ust}', 'LabelHorizontalAlignment', 'left');
xline(czas_skok, 'r:', 'Chwila skoku t_0 = 500 s', 'LabelVerticalAlignment', 'bottom');

plot([t_p1, t_p1], [y0, y_p1], 'k:', 'LineWidth', 1.2);
plot([0, t_p1], [y_p1, y_p1], 'k:', 'LineWidth', 1.2);
plot([t_p2, t_p2], [y0, y_p2], 'k:', 'LineWidth', 1.2);
plot([0, t_p2], [y_p2, y_p2], 'k:', 'LineWidth', 1.2);

h_p1 = plot(t_p1, y_p1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
h_p2 = plot(t_p2, y_p2, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');

text(t_p1+35, y_p1-0.05, sprintf('28.3%% (t = %.1f s, T_{wew} = %.2f °C)', t_p1, y_p1), ...
     'FontWeight', 'bold', 'Color', [0 0.5 0], 'FontSize', 10);
text(t_p2+35, y_p2-0.05, sprintf('63.2%% (t = %.1f s, T_{wew} = %.2f °C)', t_p2, y_p2), ...
     'FontWeight', 'bold', 'Color', [0.6 0.5 0], 'FontSize', 10);

ylim([y0 - 0.1, y_ust + 0.15]);
legend([h_ob, h_p1, h_p2], {'Odpowiedź obiektu T_{wew}', 'Punkt 28.3% (t_1 = 238.0 s)', 'Punkt 63.2% (t_2 = 522.0 s)'}, ...
       'Location', 'southeast', 'FontSize', 10);

exportgraphics(f1, 'metoda_dwupunktowa.png', 'Resolution', 300);

% --- 2. Wykres: metoda_styczna.png ---
f2 = figure('Color', 'w', 'Position', [100, 100, 950, 480]); 
hold on; grid on; box on;
title('Identyfikacja metodą stycznej - konstrukcja graficzna', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11); 
ylabel('Temperatura wewnętrzna T_{wew} [°C]', 'FontSize', 11);

h_ob2 = plot(t, aTwew, 'b-', 'LineWidth', 2.2);
yline(y0, 'k--', 'Poziom początkowy T_{wew0}', 'LabelHorizontalAlignment', 'left');
yline(y_ust, 'k--', 'Wartość ustalona T_{wew,ust}', 'LabelHorizontalAlignment', 'left');
xline(czas_skok, 'k:', 'Chwila skoku t_0 = 500 s', 'LabelVerticalAlignment', 'top');

h_st = plot([t_start_st, t_koniec_st], [y0, y_ust], 'r-', 'LineWidth', 1.8);
h_przeg = plot(t_przeg, y_przeg, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(t_start_st, y0, 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
plot(t_koniec_st, y_ust, 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

plot([t_koniec_st, t_koniec_st], [y0, y_ust], 'k:', 'LineWidth', 1.2);
plot([czas_skok, t_start_st], [y0, y0], 'm-', 'LineWidth', 3.5);
plot([t_start_st, t_koniec_st], [y0, y0], 'g-', 'LineWidth', 3.5);

text((czas_skok + t_start_st)/2, y0 + 0.08, sprintf('T_0 = %.2f s', T0_st), ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'm', 'FontSize', 10);
text((t_start_st + t_koniec_st)/2, y0 + 0.08, sprintf('T = %.2f s', T_st), ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', [0 0.5 0], 'FontSize', 10);
text(t_przeg+40, y_przeg-0.05, sprintf('Pkt przegięcia (t = %.1f s, T_{wew} = %.2f °C)', t_przeg, y_przeg), ...
     'FontWeight', 'bold', 'Color', 'r', 'FontSize', 10);

ylim([y0 - 0.12, y_ust + 0.15]);
legend([h_ob2, h_st, h_przeg], {'Odpowiedź obiektu T_{wew}', 'Styczna w pkt. przegięcia', 'Punkt przegięcia'}, ...
       'Location', 'southeast', 'FontSize', 10);

exportgraphics(f2, 'metoda_styczna.png', 'Resolution', 300);

% --- 3. Wykres: porownanie_identyfikacji.png ---
f3 = figure('Color', 'w', 'Position', [100, 100, 950, 480]); 
hold on; grid on; box on;
title('Weryfikacja modeli FOPDT - porównanie z obiektem rzeczywistym', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11); 
ylabel('Temperatura wewnętrzna T_{wew} [°C]', 'FontSize', 11);

p1 = plot(t, aTwew, 'b-', 'LineWidth', 2.5);
p2 = plot(t, Twew_model_2p, 'm-.', 'LineWidth', 2.0);
p3 = plot(t, Twew_model_st, 'g--', 'LineWidth', 2.0);

ylim([y0 - 0.05, y_ust + 0.1]);
legend([p1, p2, p3], ...
       {'Obiekt rzeczywisty (T_{wew})', ...
        sprintf('Model dwupunktowy (K=%.2f, T=%.2f s, T_0=%.2f s)', K_2p, T_2p, T0_2p), ...
        sprintf('Model ze stycznej (K=%.2f, T=%.2f s, T_0=%.2f s)', K_st, T_st, T0_st)}, ...
       'Location', 'southeast', 'FontSize', 10);

exportgraphics(f3, 'porownanie_identyfikacji.png', 'Resolution', 300);

% --- 4. Wykres: metoda_momentow.png ---
f4 = figure('Color', 'w', 'Position', [100, 100, 950, 480]); 
hold on; grid on; box on;
title('Porównanie modelu z metody momentów (II rzędu) z obiektem i modelami FOPDT', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11); 
ylabel('Temperatura wewnętrzna T_{wew} [°C]', 'FontSize', 11);

m_p1 = plot(t, aTwew, 'b-', 'LineWidth', 2.5);
m_p2 = plot(t, Twew_model_mom, 'r:', 'LineWidth', 2.2);
m_p3 = plot(t, Twew_model_2p, 'm-.', 'LineWidth', 1.5);
m_p4 = plot(t, Twew_model_st, 'g--', 'LineWidth', 1.5);

ylim([y0 - 0.05, y_ust + 0.1]);
legend([m_p1, m_p2, m_p3, m_p4], ...
       {'Obiekt rzeczywisty (T_{wew})', 'Model z metody momentów (II rzędu)', ...
        'Model dwupunktowy (FOPDT)', 'Model ze stycznej (FOPDT)'}, ...
       'Location', 'southeast', 'FontSize', 10);

exportgraphics(f4, 'metoda_momentow.png', 'Resolution', 300);

fprintf('Wszystkie 4 wykresy zostaly wyeksportowane w rozdzielczosci 300 DPI.\n');