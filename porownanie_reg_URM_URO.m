% =========================================================================
% PORÓWNANIE URM (MODEL) vs URO (OBIEKT) - SKOK +1 °C
% 4 pliki PNG (Styczna PV/CV, Dwupunktowa PV/CV), po 3 panele (3x1) na plik
% Każdy panel = 1 metoda strojenia (Z-N, CHR 20%, pidtune)
% =========================================================================
clear all; close all; clc;

% 1. INICJALIZACJA OBIEKTU FIZYCZNEGO (dla URO.slx)
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

% Parametry eksperymentu (skok +1 °C)
czas      = 4500;
czas_skok = 500;
dSPc_reg  = 1.0;

% Zerowanie zakłóceń i skoku w KPTz
dSPc = 0; dTzew = 0; dFg = 0; dFw = 0; dQt = 0;
wp = Fw0;

% 2. MODELE Z IDENTYFIKACJI I NASTAWY PI
K_st = 155.42; T_st = 618.08; T0_st = 62.15;
K_2p = 155.42; T_2p = 426.00; T0_2p = 96.00;

s = tf('s');
G_st = (K_st / (T_st * s + 1)) * exp(-T0_st * s);
G_2p = (K_2p / (T_2p * s + 1)) * exp(-T0_2p * s);

% Nastawy ze stycznej
Kp_st_ZN  = (0.9 * T_st) / (K_st * T0_st);  Ti_st_ZN  = 3.33 * T0_st;
Kp_st_CHR = (0.6 * T_st) / (K_st * T0_st);  Ti_st_CHR = 1.0  * T_st;
C_st_pt   = pidtune(G_st, 'PI');
Kp_st_pt  = C_st_pt.Kp;                     Ti_st_pt  = C_st_pt.Kp / C_st_pt.Ki;

% Nastawy z modelu dwupunktowego
Kp_2p_ZN  = (0.9 * T_2p) / (K_2p * T0_2p);  Ti_2p_ZN  = 3.33 * T0_2p;
Kp_2p_CHR = (0.6 * T_2p) / (K_2p * T0_2p);  Ti_2p_CHR = 1.0  * T_2p;
C_2p_pt   = pidtune(G_2p, 'PI');
Kp_2p_pt  = C_2p_pt.Kp;                     Ti_2p_pt  = C_2p_pt.Kp / C_2p_pt.Ki;

% {Model, Metoda, Nr_panelu (1..3), Kp, Ti, Transmitancja_URM}
konf = {
    'Styczna',     'Ziegler-Nichols', 1, Kp_st_ZN,  Ti_st_ZN,  G_st;
    'Styczna',     'CHR (20%)',       2, Kp_st_CHR, Ti_st_CHR, G_st;
    'Styczna',     'pidtune',         3, Kp_st_pt,  Ti_st_pt,  G_st;
    'Dwupunktowy', 'Ziegler-Nichols', 1, Kp_2p_ZN,  Ti_2p_ZN,  G_2p;
    'Dwupunktowy', 'CHR (20%)',       2, Kp_2p_CHR, Ti_2p_CHR, G_2p;
    'Dwupunktowy', 'pidtune',         3, Kp_2p_pt,  Ti_2p_pt,  G_2p
};

% 3. UTWORZENIE 4 FIGUR Z UKŁADEM 3x1 (PO 3 METODY NA JEDEN PNG)
f_st_PV = figure('Name', 'Styczna - PV (URM vs URO)', 'Color', 'w', 'Position', [40,  40, 860, 820]);
tl_st_PV = tiledlayout(f_st_PV, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_st_PV, 'Porównanie URM vs URO (nastawy ze stycznej) – wielkość regulowana PV (T_{wew})', 'FontSize', 12, 'FontWeight', 'bold');

f_st_CV = figure('Name', 'Styczna - CV (URM vs URO)', 'Color', 'w', 'Position', [920, 40, 860, 820]);
tl_st_CV = tiledlayout(f_st_CV, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_st_CV, 'Porównanie URM vs URO (nastawy ze stycznej) – sygnał sterujący CV (F_w)', 'FontSize', 12, 'FontWeight', 'bold');

f_2p_PV = figure('Name', 'Dwupunktowa - PV (URM vs URO)', 'Color', 'w', 'Position', [80,  40, 860, 820]);
tl_2p_PV = tiledlayout(f_2p_PV, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_2p_PV, 'Porównanie URM vs URO (nastawy z m. dwupunktowej) – wielkość regulowana PV (T_{wew})', 'FontSize', 12, 'FontWeight', 'bold');

f_2p_CV = figure('Name', 'Dwupunktowa - CV (URM vs URO)', 'Color', 'w', 'Position', [960, 40, 860, 820]);
tl_2p_CV = tiledlayout(f_2p_CV, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_2p_CV, 'Porównanie URM vs URO (nastawy z m. dwupunktowej) – sygnał sterujący CV (F_w)', 'FontSize', 12, 'FontWeight', 'bold');

% Kolory: URO (obiekt) = niebieski ciągły, URM (model) = czerwony przerywany
kol_URO = [0.00, 0.45, 0.74];
kol_URM = [0.85, 0.20, 0.15];

t_urm = (0:1:czas)';
idx_s = t_urm >= czas_skok;
t_sp  = [0, czas_skok, czas_skok, czas];
y_sp  = [Twew0, Twew0, Twew0 + dSPc_reg, Twew0 + dSPc_reg];

n = size(konf, 1);
W = struct();

% 4. PĘTLA SYMULACJI I RYSOWANIA NA OSOBNYCH PANELACH
for i = 1:n
    mod_n = konf{i, 1};
    met_n = konf{i, 2};
    nr_p  = konf{i, 3};
    Kp    = konf{i, 4};
    Ti    = konf{i, 5};
    Gm    = konf{i, 6};
    
    % --- A. URM (model liniowy FOPDT) ---
    R_PI = Kp * (1 + 1 / (Ti * s));
    PV_u = Twew0 * ones(size(t_urm));
    CV_u = Fw0   * ones(size(t_urm));
    PV_u(idx_s) = Twew0 + dSPc_reg * step(feedback(R_PI * Gm, 1), t_urm(idx_s) - czas_skok);
    CV_u(idx_s) = Fw0   + dSPc_reg * step(feedback(R_PI, Gm),     t_urm(idx_s) - czas_skok);
    [tr_u, ise_u, iae_u] = licz_wskazniki(t_urm, PV_u, czas_skok, Twew0 + dSPc_reg, dSPc_reg, 0.05);
    
    % --- B. URO (obiekt fizyczny URO.slx) ---
    simOut = sim('URO', 'StopTime', num2str(czas), ...
                 'ReturnWorkspaceOutputs', 'on', 'SaveTime', 'on', 'TimeSaveName', 'tout');
    [t_o, PV_o] = wyciagnij(simOut, 'aTwew');
    [~,   CV_o] = wyciagnij(simOut, 'aFw');
    [tr_o, ise_o, iae_o] = licz_wskazniki(t_o, PV_o, czas_skok, Twew0 + dSPc_reg, dSPc_reg, 0.05);
    
    % Zapis wyników do tabeli
    W(i).model = mod_n; W(i).metoda = met_n; W(i).Kp = Kp; W(i).Ti = Ti;
    W(i).tr_u = tr_u; W(i).ise_u = ise_u; W(i).iae_u = iae_u;
    W(i).tr_o = tr_o; W(i).ise_o = ise_o; W(i).iae_o = iae_o;
    
    % Wybór odpowiedniego układu tiledlayout
    if strcmp(mod_n, 'Styczna')
        tlPV = tl_st_PV; tlCV = tl_st_CV;
    else
        tlPV = tl_2p_PV; tlCV = tl_2p_CV;
    end
    
    tytul_panelu = sprintf('%s (K_p = %.5f, T_i = %.2f s)', met_n, Kp, Ti);
    
    % --- PANEL PV ---
    axPV = nexttile(tlPV, nr_p); hold(axPV, 'on'); grid(axPV, 'on'); box(axPV, 'on');
    plot(axPV, t_o,   PV_o, '-',  'Color', kol_URO, 'LineWidth', 2.0, 'DisplayName', 'URO (obiekt rzeczywisty)');
    plot(axPV, t_urm, PV_u, '--', 'Color', kol_URM, 'LineWidth', 1.8, 'DisplayName', 'URM (model liniowy)');
    plot(axPV, t_sp,  y_sp, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Wartość zadana SP');
    yline(axPV, Twew0 + 1.05*dSPc_reg, 'k--', 'Alpha', 0.3, 'HandleVisibility', 'off');
    yline(axPV, Twew0 + 0.95*dSPc_reg, 'k--', 'Alpha', 0.3, 'HandleVisibility', 'off');
    title(axPV, tytul_panelu, 'FontSize', 10.5);
    ylabel(axPV, 'T_{wew} [°C]');
    xlim(axPV, [400, 3500]);
    if nr_p == 3, xlabel(axPV, 'Czas t [s]'); end
    legend(axPV, 'show', 'Location', 'southeast', 'FontSize', 8.5);
    
    % --- PANEL CV ---
    axCV = nexttile(tlCV, nr_p); hold(axCV, 'on'); grid(axCV, 'on'); box(axCV, 'on');
    plot(axCV, t_o,   CV_o, '-',  'Color', kol_URO, 'LineWidth', 2.0, 'DisplayName', 'URO (obiekt rzeczywisty)');
    plot(axCV, t_urm, CV_u, '--', 'Color', kol_URM, 'LineWidth', 1.8, 'DisplayName', 'URM (model liniowy)');
    yline(axCV, Fw0, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Punkt pracy F_{w0}');
    title(axCV, tytul_panelu, 'FontSize', 10.5);
    ylabel(axCV, 'F_w [kg/s]');
    xlim(axCV, [400, 3500]);
    if nr_p == 3, xlabel(axCV, 'Czas t [s]'); end
    legend(axCV, 'show', 'Location', 'northeast', 'FontSize', 8.5);
end

% 5. EKSPORT 4 PLIKÓW PNG (300 DPI)
exportgraphics(f_st_PV, 'Porownanie_styczna_PV.png',     'Resolution', 300);
exportgraphics(f_st_CV, 'Porownanie_styczna_CV.png',     'Resolution', 300);
exportgraphics(f_2p_PV, 'Porownanie_dwupunktowa_PV.png', 'Resolution', 300);
exportgraphics(f_2p_CV, 'Porownanie_dwupunktowa_CV.png', 'Resolution', 300);

% 6. TABELA WYNIKÓW W KONSOLI I GOTOWY KOD LATEX
fprintf('\n==============================================================================================================\n');
fprintf('                      ZESTAWIENIE WYNIKÓW URM vs URO (SKOK +1 °C)\n');
fprintf('==============================================================================================================\n');
fprintf('%-11s | %-15s | %-7s | %-6s | %-8s | %-8s | %-7s | %-7s | %-7s | %-7s\n', ...
        'Model', 'Metoda', 'Kp', 'Ti[s]', 'tr URM', 'tr URO', 'ISE URM', 'ISE URO', 'IAE URM', 'IAE URO');
fprintf('--------------------------------------------------------------------------------------------------------------\n');
for k = 1:n
    fprintf('%-11s | %-15s | %7.5f | %6.1f | %8.1f | %8.1f | %7.2f | %7.2f | %7.2f | %7.2f\n', ...
            W(k).model, W(k).metoda, W(k).Kp, W(k).Ti, ...
            W(k).tr_u, W(k).tr_o, W(k).ise_u, W(k).ise_o, W(k).iae_u, W(k).iae_o);
end
fprintf('==============================================================================================================\n');

fprintf('\n--- GOTOWA TABELA ZBIORCZA DO LATEXA (URM vs URO) ---\n');
fprintf('\\begin{table}[H]\n\\centering\n\\small\n');
fprintf('\\caption{Porównanie wskaźników jakości regulacji PI dla URM i URO (skok $\\Delta SP = +1\\,^\\circ\\mathrm{C}$)}\n');
fprintf('\\label{tab:porownanie_urm_uro}\n');
fprintf('\\begin{tabular}{|l|l|c|c|c|c|c|c|c|c|}\n\\hline\n');
fprintf('\\textbf{Model} & \\textbf{Metoda} & \\(K_p\\) & \\(T_i\\,[\\mathrm{s}]\\) & \\(t_{r,\\mathrm{URM}}\\) & \\(t_{r,\\mathrm{URO}}\\) & \\(\\mathrm{ISE}_{\\mathrm{URM}}\\) & \\(\\mathrm{ISE}_{\\mathrm{URO}}\\) & \\(\\mathrm{IAE}_{\\mathrm{URM}}\\) & \\(\\mathrm{IAE}_{\\mathrm{URO}}\\) \\\\ \\hline\n');
for k = 1:n
    fprintf('%s & %s & %.5f & %.1f & %.1f & %.1f & %.2f & %.2f & %.2f & %.2f \\\\ \\hline\n', ...
            W(k).model, W(k).metoda, W(k).Kp, W(k).Ti, ...
            W(k).tr_u, W(k).tr_o, W(k).ise_u, W(k).ise_o, W(k).iae_u, W(k).iae_o);
end
fprintf('\\end{tabular}\n\\end{table}\n');

% =========================================================================
% FUNKCJE POMOCNICZE
% =========================================================================
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

function [tr, ISE, IAE] = licz_wskazniki(t, y, t_skok, SP_doc, dSP, eps_proc)
    idx = t >= t_skok;
    t_r = t(idx) - t_skok;
    e   = SP_doc - y(idx);
    
    ISE = trapz(t_r, e.^2);
    IAE = trapz(t_r, abs(e));
    
    idx_poza = find(abs(e) > eps_proc * abs(dSP), 1, 'last');
    if isempty(idx_poza)
        tr = 0;
    elseif idx_poza == length(t_r)
        tr = NaN;
    else
        tr = t_r(idx_poza);
    end
end