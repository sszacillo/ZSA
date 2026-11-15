clear all; clc;
load('modele_transmitancje.mat');

figure('Name', 'Metoda stycznej', 'Color', 'w');
hold on; grid on;
title('Identyfikacja metodą stycznej - konstrukcja graficzna');
xlabel('Czas t [s]');
ylabel('Temperatura wewnętrzna T_{wew} [°C]');

% 1. Odpowiedź obiektu i asymptoty
h_ob = plot(t, aTwew, 'b-', 'LineWidth', 2.2);
yline(y0, 'k--', 'Poziom początkowy T_{wew0}', 'LabelHorizontalAlignment', 'left');
yline(y_ust, 'k--', 'Wartość ustalona T_{wew,ust}', 'LabelHorizontalAlignment', 'left');
xline(czas_skok, 'k:', 'Chwila skoku t_0 = 500 s', 'LabelVerticalAlignment', 'top');

% 2. Linia stycznej w punkcie przegięcia
h_st = plot([t_start_st, t_koniec_st], [y0, y_ust], 'r-', 'LineWidth', 1.8);

% 3. Punkt przegięcia i punkty przecięcia z asymptotami
h_przeg = plot(t_przeg, y_przeg, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(t_start_st, y0, 'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(t_koniec_st, y_ust, 'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'k');

% 4. Rzutowanie górnego przecięcia na poziom początkowy y0
plot([t_koniec_st, t_koniec_st], [y0, y_ust], 'k:', 'LineWidth', 1.2);

% 5. Wektory/odcinki pokazujące T0 oraz T na osi poziomej y0
plot([czas_skok, t_start_st], [y0, y0], 'm-', 'LineWidth', 3);
plot([t_start_st, t_koniec_st], [y0, y0], 'g-', 'LineWidth', 3);

text((czas_skok + t_start_st)/2, y0 + 0.07, sprintf('T_0 = %.1f s', T0_st), ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'm');
text((t_start_st + t_koniec_st)/2, y0 + 0.07, sprintf('T = %.1f s', T_st), ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', [0 0.5 0]);
text(t_przeg + 40, y_przeg - 0.05, sprintf('Punkt przegięcia\n(t = %.1f s, T_{wew} = %.2f °C)', t_przeg, y_przeg), ...
     'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');

ylim([y0 - 0.12, y_ust + 0.15]);
legend([h_ob, h_st, h_przeg], ...
       {'Odpowiedź obiektu T_{wew}', 'Styczna w punkcie przegięcia', 'Punkt przegięcia'}, ...
       'FontSize', 10, 'Location', 'southeast');