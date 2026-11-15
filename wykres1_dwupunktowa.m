clear all; clc;
load('modele_transmitancje.mat');

figure('Name', 'Metoda dwupunktowa', 'Color', 'w');
hold on; grid on;
title('Identyfikacja metodą dwupunktową - konstrukcja graficzna');
xlabel('Czas t [s]');
ylabel('Temperatura wewnętrzna T_{wew} [°C]');

% 1. Odpowiedź obiektu i poziomy odniesienia
h_ob = plot(t, aTwew, 'b-', 'LineWidth', 2.2);
yline(y0, 'k--', 'Poziom początkowy T_{wew0}', 'LabelHorizontalAlignment', 'left');
yline(y_ust, 'k--', 'Wartość ustalona T_{wew,ust}', 'LabelHorizontalAlignment', 'left');
xline(czas_skok, 'r:', 'Chwila skoku t_0 = 500 s', 'LabelVerticalAlignment', 'bottom');

% 2. Linie rzutujące punkty na osie
plot([t_p1, t_p1], [y0, y_p1], 'k:', 'LineWidth', 1.2);
plot([0, t_p1], [y_p1, y_p1], 'k:', 'LineWidth', 1.2);
plot([t_p2, t_p2], [y0, y_p2], 'k:', 'LineWidth', 1.2);
plot([0, t_p2], [y_p2, y_p2], 'k:', 'LineWidth', 1.2);

% 3. Punkty pomiarowe 28.3% i 63.2%
h_p1 = plot(t_p1, y_p1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
h_p2 = plot(t_p2, y_p2, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');

% 4. Podpisy na wykresie
text(t_p1 + 35, y_p1 - 0.04, sprintf('28.3%% (t = %.1f s, T_{wew} = %.2f °C)', t_p1, y_p1), ...
     'FontSize', 10, 'FontWeight', 'bold', 'Color', [0 0.5 0]);
text(t_p2 + 35, y_p2 - 0.04, sprintf('63.2%% (t = %.1f s, T_{wew} = %.2f °C)', t_p2, y_p2), ...
     'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.6 0.5 0]);

ylim([y0 - 0.1, y_ust + 0.15]);
legend([h_ob, h_p1, h_p2], ...
       {'Odpowiedź obiektu T_{wew}', 'Punkt 28.3% przyrostu', 'Punkt 63.2% przyrostu'}, ...
       'FontSize', 10, 'Location', 'southeast');