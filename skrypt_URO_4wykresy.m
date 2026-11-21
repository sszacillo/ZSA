% =========================================================================
% SYMULACJA URO (URO.slx) - 4 OSOBNE WYKRESY (PV i CV: Styczna / 2-pkt)
% Metody doboru nastaw: Ziegler-Nichols, CHR (20%), pidtune
% =========================================================================
clear all; close all; clc;

% 1. INICJALIZACJA OBIEKTU FIZYCZNEGO (Wezel1Pom1CvwCvg / URO.slx)
aparam1;
QgN = 2000; TzewN = -20; TwewN = 20;
TgzN = 90;  TgpN = 70;   TwzN = 130; TwpN = 110;

FgN   = QgN / (cpw * (TgzN - TgpN));
Kcg   = QgN / (TgpN - TwewN);
Kcw   = QgN / (TwewN - TzewN);
FwN   = QgN / (TwzN - TwpN) / cpw;
Kco   = QgN / (TwpN - TgzN);
atz_c = (TwzN - TzewN) / (TwewN - TzewN);
btz_c = (TwzN - TwewN) / (TwewN - TzewN);

SPc0  = TwewN; Tzew0 = TzewN;
Twz0  = atz_c * SPc0 - btz_c * Tzew0;
Fg0   = FgN;   Fw0   = FwN;   Qt0 = 0;

A = [(-Kcg-Kcw)/Cvws,  Kcg/Cvws,                 0,                  0;
      Kcg/Cvg,         (-cpw*Fg0-Kcg)/Cvg,       cpw*Fg0/Cvg,        0;
      0,               cpw*Fg0/Cvo,              (-cpw*Fg0-Kco)/Cvo, Kco/Cvo;
      0,               0,                        Kco/Cvo,            (-cpw*Fw0-Kco)/Cvo];
B = [0, Kcw/Cvws; 0, 0; 0, 0; cpw*Fw0/Cvo, 0];
x0 = -A^(-1) * B * [Twz0; Tzew0];
Twew0 = x0(1); Tgp0 = x0(2); Tgz0 = x0(3); Twp0 = x0(4);

% Parametry eksperymentu URO
czas      = 4500;   % Czas symulacji [s]
czas_skok = 500;    % Chwila skoku wartości zadanej [s]
dSPc_reg  = 1.0;    % Skok wartości zadanej regulatora (+1 °C: 20 -> 21 °C)

% Zerowanie zakłóceń i skoku dla KPTz (aby Twz pozostało stałe)
dSPc = 0; dTzew = 0; dFg = 0; dFw = 0; dQt = 0;

% Warunek początkowy integratora 1/s w PI regulator (bezuderzeniowy start)
wp = Fw0;

% 2. PARAMETRY MODELI Z IDENTYFIKACJI I WYZNACZENIE NASTAW PI
K_st = 155.42; T_st = 618.08; T0_st = 62.15;
K_2p = 155.42; T_2p = 426.00; T0_2p = 96.00;

s = tf('s');
G_st = (K_st / (T_st * s + 1)) * exp(-T0_st * s);
G_2p = (K_2p / (T_2p * s + 1)) * exp(-T0_2p * s);

% Nastawy z modelu ze stycznej
Kp_st_ZN  = (0.9 * T_st) / (K_st * T0_st);  Ti_st_ZN  = 3.33 * T0_st;
Kp_st_CHR = (0.6 * T_st) / (K_st * T0_st);  Ti_st_CHR = 1.0  * T_st;
C_st_pt   = pidtune(G_st, 'PI');
Kp_st_pt  = C_st_pt.Kp;                     Ti_st_pt  = C_st_pt.Kp / C_st_pt.Ki;

% Nastawy z modelu dwupunktowego
Kp_2p_ZN  = (0.9 * T_2p) / (K_2p * T0_2p);  Ti_2p_ZN  = 3.33 * T0_2p;
Kp_2p_CHR = (0.6 * T_2p) / (K_2p * T0_2p);  Ti_2p_CHR = 1.0  * T_2p;
C_2p_pt   = pidtune(G_2p, 'PI');
Kp_2p_pt  = C_2p_pt.Kp;                     Ti_2p_pt  = C_2p_pt.Kp / C_2p_pt.Ki;

konf = {
    'Styczna',     'Ziegler-Nichols', Kp_st_ZN,  Ti_st_ZN,  '-',  [0.85 0.15 0.15];
    'Styczna',     'CHR (20%)',       Kp_st_CHR, Ti_st_CHR, '--', [0.00 0.45 0.74];
    'Styczna',     'pidtune',         Kp_st_pt,  Ti_st_pt,  '-.', [0.15 0.65 0.20];
    'Dwupunktowy', 'Ziegler-Nichols', Kp_2p_ZN,  Ti_2p_ZN,  '-',  [0.85 0.15 0.15];
    'Dwupunktowy', 'CHR (20%)',       Kp_2p_CHR, Ti_2p_CHR, '--', [0.00 0.45 0.74];
    'Dwupunktowy', 'pidtune',         Kp_2p_pt,  Ti_2p_pt,  '-.', [0.15 0.65 0.20]
};

% 3. PRZYGOTOWANIE 4 OSOBNYCH WYKRESÓW DLA URO
f_st_PV = figure('Name', 'URO Styczna - PV', 'Color', 'w', 'Position', [50,  420, 880, 440]);
hold on; grid on; box on;
title('URO (obiekt, nastawy ze stycznej) - wielkość regulowana PV (T_{wew})', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11); ylabel('Temperatura T_{wew} [°C]', 'FontSize', 11);

f_st_CV = figure('Name', 'URO Styczna - CV', 'Color', 'w', 'Position', [960, 420, 880, 440]);
hold on; grid on; box on;
title('URO (obiekt, nastawy ze stycznej) - sygnał sterujący CV (F_w)', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11); ylabel('Przepływ F_w [kg/s]', 'FontSize', 11);

f_2p_PV = figure('Name', 'URO Dwupunktowa - PV', 'Color', 'w', 'Position', [50,  40, 880, 440]);
hold on; grid on; box on;
title('URO (obiekt, nastawy z m. dwupunktowej) - wielkość regulowana PV (T_{wew})', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11); ylabel('Temperatura T_{wew} [°C]', 'FontSize', 11);

f_2p_CV = figure('Name', 'URO Dwupunktowa - CV', 'Color', 'w', 'Position', [960, 40, 880, 440]);
hold on; grid on; box on;
title('URO (obiekt, nastawy z m. dwupunktowej) - sygnał sterujący CV (F_w)', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11); ylabel('Przepływ F_w [kg/s]', 'FontSize', 11);

% 4. PĘTLA SYMULACJI NA PLIKU "URO.slx"
n = size(konf, 1);
wyn_URO = struct();

for i = 1:n
    mod_n = konf{i, 1};
    met_n = konf{i, 2};
    Kp    = konf{i, 3};
    Ti    = konf{i, 4};
    styl  = konf{i, 5};
    kolor = konf{i, 6};
    opis  = sprintf('%s (K_p = %.5f, T_i = %.2f s)', met_n, Kp, Ti);
    
    simOut = sim('URO', 'StopTime', num2str(czas), ...
                 'ReturnWorkspaceOutputs', 'on', 'SaveTime', 'on', 'TimeSaveName', 'tout');
             
    [t_o, PV_o] = wyciagnij(simOut, 'aTwew');
    [~,   CV_o] = wyciagnij(simOut, 'aFw');
    
    % Obliczenie wskaźników jakości (tr +/- 5%, ISE, IAE)
    idx_r = t_o >= czas_skok;
    t_r   = t_o(idx_r) - czas_skok;
    e_r   = (Twew0 + dSPc_reg) - PV_o(idx_r);
    
    ISE = trapz(t_r, e_r.^2);
    IAE = trapz(t_r, abs(e_r));
    idx_poza = find(abs(e_r) > 0.05 * abs(dSPc_reg), 1, 'last');
    if isempty(idx_poza), tr = 0; elseif idx_poza == length(t_r), tr = NaN; else, tr = t_r(idx_poza); end
    
    wyn_URO(i).model = mod_n; wyn_URO(i).metoda = met_n;
    wyn_URO(i).Kp = Kp; wyn_URO(i).Ti = Ti;
    wyn_URO(i).tr = tr; wyn_URO(i).ISE = ISE; wyn_URO(i).IAE = IAE;
    
    if strcmp(mod_n, 'Styczna')
        figure(f_st_PV); plot(t_o, PV_o, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.2, 'DisplayName', opis);
        figure(f_st_CV); plot(t_o, CV_o, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.2, 'DisplayName', opis);
    else
        figure(f_2p_PV); plot(t_o, PV_o, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.2, 'DisplayName', opis);
        figure(f_2p_CV); plot(t_o, CV_o, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.2, 'DisplayName', opis);
    end
end

% 5. WYKOŃCZENIE I EKSPORT 4 WYKRESÓW (300 DPI)
t_sp = [0, czas_skok, czas_skok, czas];
y_sp = [Twew0, Twew0, Twew0 + dSPc_reg, Twew0 + dSPc_reg];

figure(f_st_PV);
plot(t_sp, y_sp, 'k:', 'LineWidth', 1.8, 'DisplayName', 'Wartość zadana SP');
xlim([400, czas]); legend('show', 'Location', 'southeast', 'FontSize', 9.5);
exportgraphics(f_st_PV, 'URO_styczna_PV.png', 'Resolution', 300);

figure(f_st_CV);
yline(Fw0, 'k:', 'Punkt pracy F_{w0}', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlim([400, czas]); legend('show', 'Location', 'northeast', 'FontSize', 9.5);
exportgraphics(f_st_CV, 'URO_styczna_CV.png', 'Resolution', 300);

figure(f_2p_PV);
plot(t_sp, y_sp, 'k:', 'LineWidth', 1.8, 'DisplayName', 'Wartość zadana SP');
xlim([400, czas]); legend('show', 'Location', 'southeast', 'FontSize', 9.5);
exportgraphics(f_2p_PV, 'URO_dwupunktowa_PV.png', 'Resolution', 300);

figure(f_2p_CV);
yline(Fw0, 'k:', 'Punkt pracy F_{w0}', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlim([400, czas]); legend('show', 'Location', 'northeast', 'FontSize', 9.5);
exportgraphics(f_2p_CV, 'URO_dwupunktowa_CV.png', 'Resolution', 300);

% 6. TABELA WYNIKÓW I KOD LATEX
fprintf('\n=================================================================================\n');
fprintf('                      WYNIKI REGULACJI: URO (Obiekt rzeczywisty)\n');
fprintf('=================================================================================\n');
fprintf('%-12s | %-16s | %-9s | %-8s | %-10s | %-8s | %-8s\n', ...
        'Model', 'Metoda', 'Kp', 'Ti [s]', 'tr(5%)[s]', 'ISE', 'IAE');
fprintf('---------------------------------------------------------------------------------\n');
for k = 1:n
    fprintf('%-12s | %-16s | %9.5f | %8.2f | %10.1f | %8.2f | %8.2f\n', ...
            wyn_URO(k).model, wyn_URO(k).metoda, wyn_URO(k).Kp, wyn_URO(k).Ti, ...
            wyn_URO(k).tr, wyn_URO(k).ISE, wyn_URO(k).IAE);
end
fprintf('=================================================================================\n');

fprintf('\n--- GOTOWY KOD TABELI DO LATEXA (URO) ---\n');
fprintf('\\begin{table}[H]\n\\centering\n');
fprintf('\\caption{Zestawienie nastaw regulatora PI i wskaźników jakości dla URO}\n');
fprintf('\\label{tab:uro_wyniki}\n');
fprintf('\\begin{tabular}{|l|l|c|c|c|c|c|}\n\\hline\n');
fprintf('\\textbf{Model} & \\textbf{Metoda} & \\(K_p\\) & \\(T_i\\,[\\mathrm{s}]\\) & \\(t_r\\,[\\mathrm{s}]\\) & \\textbf{ISE} & \\textbf{IAE} \\\\ \\hline\n');
for k = 1:n
    fprintf('%s & %s & %.5f & %.2f & %.1f & %.2f & %.2f \\\\ \\hline\n', ...
            wyn_URO(k).model, wyn_URO(k).metoda, wyn_URO(k).Kp, wyn_URO(k).Ti, ...
            wyn_URO(k).tr, wyn_URO(k).ISE, wyn_URO(k).IAE);
end
fprintf('\\end{tabular}\n\\end{table}\n');

function [t_out, y_out] = wyciagnij(simOut, nazwa)
    dane = simOut.get(nazwa);
    if isa(dane, 'timeseries')
        t_out = dane.Time; y_out = squeeze(dane.Data);
    elseif isstruct(dane) && isfield(dane, 'signals')
        t_out = dane.time; y_out = squeeze(dane.signals.values);
    else
        t_out = simOut.tout; y_out = squeeze(dane);
    end
end