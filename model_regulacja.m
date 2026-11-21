% =========================================================================
% SYMULACJA URM (model.mdl) - 4 OSOBNE WYKRESY (PV i CV: Styczna / 2-pkt)
% Metody doboru nastaw: Ziegler-Nichols, CHR (20%), pidtune
% =========================================================================
clear all; close all; clc;

% 1. Parametry modeli z identyfikacji
K_st_baza = 155.42; T_st_baza = 618.08; T0_st_baza = 62.15;
K_2p_baza = 155.42; T_2p_baza = 426.00; T0_2p_baza = 96.00;

% Punkt pracy i parametry eksperymentu
Twew0     = 20.0;       % Temperatura początkowa [°C]
Fw0       = 0.023866;   % Nominalny przepływ FwN [kg/s]
czas      = 4500;       % Czas symulacji [s]
czas_skok = 500;        % Chwila skoku wartości zadanej [s]
dSPc      = 1.0;        % Skok wartości zadanej (+1 °C: 20 -> 21 °C)

% 2. Wyznaczenie nastaw regulatora PI (Ziegler-Nichols, CHR 20%, pidtune)
s = tf('s');
G_st = (K_st_baza / (T_st_baza * s + 1)) * exp(-T0_st_baza * s);
G_2p = (K_2p_baza / (T_2p_baza * s + 1)) * exp(-T0_2p_baza * s);

% Nastawy na bazie modelu ze stycznej
Kp_st_ZN  = (0.9 * T_st_baza) / (K_st_baza * T0_st_baza);  Ti_st_ZN  = 3.33 * T0_st_baza;
Kp_st_CHR = (0.6 * T_st_baza) / (K_st_baza * T0_st_baza);  Ti_st_CHR = 1.0  * T_st_baza;
C_st_pt   = pidtune(G_st, 'PI');
Kp_st_pt  = C_st_pt.Kp;                                    Ti_st_pt  = C_st_pt.Kp / C_st_pt.Ki;

% Nastawy na bazie modelu dwupunktowego
Kp_2p_ZN  = (0.9 * T_2p_baza) / (K_2p_baza * T0_2p_baza);  Ti_2p_ZN  = 3.33 * T0_2p_baza;
Kp_2p_CHR = (0.6 * T_2p_baza) / (K_2p_baza * T0_2p_baza);  Ti_2p_CHR = 1.0  * T_2p_baza;
C_2p_pt   = pidtune(G_2p, 'PI');
Kp_2p_pt  = C_2p_pt.Kp;                                    Ti_2p_pt  = C_2p_pt.Kp / C_2p_pt.Ki;

% Konfiguracja 6 symulacji: {Model, Metoda, Kp, Ti, K_ob, T_ob, T0_ob, Styl, Kolor}
konf = {
    'Styczna',     'Ziegler-Nichols', Kp_st_ZN,  Ti_st_ZN,  K_st_baza, T_st_baza, T0_st_baza, '-',  [0.85 0.15 0.15];
    'Styczna',     'CHR (20%)',       Kp_st_CHR, Ti_st_CHR, K_st_baza, T_st_baza, T0_st_baza, '--', [0.00 0.45 0.74];
    'Styczna',     'pidtune',         Kp_st_pt,  Ti_st_pt,  K_st_baza, T_st_baza, T0_st_baza, '-.', [0.15 0.65 0.20];
    'Dwupunktowy', 'Ziegler-Nichols', Kp_2p_ZN,  Ti_2p_ZN,  K_2p_baza, T_2p_baza, T0_2p_baza, '-',  [0.85 0.15 0.15];
    'Dwupunktowy', 'CHR (20%)',       Kp_2p_CHR, Ti_2p_CHR, K_2p_baza, T_2p_baza, T0_2p_baza, '--', [0.00 0.45 0.74];
    'Dwupunktowy', 'pidtune',         Kp_2p_pt,  Ti_2p_pt,  K_2p_baza, T_2p_baza, T0_2p_baza, '-.', [0.15 0.65 0.20]
};

% 3. Przygotowanie 4 osobnych okien wykresów
f_st_PV = figure('Name', 'Styczna - PV', 'Color', 'w', 'Position', [50,  420, 880, 440]);
hold on; grid on; box on;
title('URM (model ze stycznej) - wielkość regulowana PV (T_{wew})', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11); ylabel('Temperatura T_{wew} [°C]', 'FontSize', 11);

f_st_CV = figure('Name', 'Styczna - CV', 'Color', 'w', 'Position', [960, 420, 880, 440]);
hold on; grid on; box on;
title('URM (model ze stycznej) - sygnał sterujący CV (F_w)', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11); ylabel('Przepływ F_w [kg/s]', 'FontSize', 11);

f_2p_PV = figure('Name', 'Dwupunktowa - PV', 'Color', 'w', 'Position', [50,  40, 880, 440]);
hold on; grid on; box on;
title('URM (model dwupunktowy) - wielkość regulowana PV (T_{wew})', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11); ylabel('Temperatura T_{wew} [°C]', 'FontSize', 11);

f_2p_CV = figure('Name', 'Dwupunktowa - CV', 'Color', 'w', 'Position', [960, 40, 880, 440]);
hold on; grid on; box on;
title('URM (model dwupunktowy) - sygnał sterujący CV (F_w)', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11); ylabel('Przepływ F_w [kg/s]', 'FontSize', 11);

% 4. Pętla symulacji (z automatyczną obsługą model.mdl oraz zabezpieczeniem)
n = size(konf, 1);
wyn_URM = struct();
uzyto_simulinka = false;

for i = 1:n
    mod_n = konf{i, 1};
    met_n = konf{i, 2};
    Kp    = konf{i, 3};
    Ti    = konf{i, 4};
    
    % Zmienne używane w blokach Transfer Fcn i Transport Delay2 w model.mdl
    K_2p  = konf{i, 5};
    T_2p  = konf{i, 6};
    T0_2p = konf{i, 7};
    
    styl  = konf{i, 8};
    kolor = konf{i, 9};
    opis  = sprintf('%s (K_p = %.5f, T_i = %.2f s)', met_n, Kp, Ti);
    
    % Próba symulacji bezpośrednio w model.mdl (wymaga bloków To Workspace: PV_urm, dCV_urm)
    udalo_sie_sim = false;
    if exist('model.mdl', 'file') || exist('model.slx', 'file')
        try
            simOut = sim('model', 'StopTime', num2str(czas), ...
                         'ReturnWorkspaceOutputs', 'on', 'SaveTime', 'on', 'TimeSaveName', 'tout');
            [t_u, PV_u]  = wyciagnij(simOut, 'PV_urm');
            [~,   dCV_u] = wyciagnij(simOut, 'dCV_urm');
            CV_u = Fw0 + dCV_u;
            udalo_sie_sim = true;
            uzyto_simulinka = true;
        catch
            udalo_sie_sim = false;
        end
    end
    
    % Równoważne wyliczenie z transmitancji, jeśli brak bloków To Workspace w model.mdl
    if ~udalo_sie_sim
        t_u   = (0:1:czas)';
        idx_s = t_u >= czas_skok;
        Gm    = (K_2p / (T_2p * s + 1)) * exp(-T0_2p * s);
        R_PI  = Kp * (1 + 1 / (Ti * s));
        
        PV_u = Twew0 * ones(size(t_u));
        CV_u = Fw0   * ones(size(t_u));
        PV_u(idx_s) = Twew0 + dSPc * step(feedback(R_PI * Gm, 1), t_u(idx_s) - czas_skok);
        CV_u(idx_s) = Fw0   + dSPc * step(feedback(R_PI, Gm),     t_u(idx_s) - czas_skok);
    end
    
    % Obliczenie wskaźników jakości (tr dla +/- 5%, ISE, IAE)
    idx_r = t_u >= czas_skok;
    t_r   = t_u(idx_r) - czas_skok;
    e_r   = (Twew0 + dSPc) - PV_u(idx_r);
    
    ISE = trapz(t_r, e_r.^2);
    IAE = trapz(t_r, abs(e_r));
    idx_poza = find(abs(e_r) > 0.05 * abs(dSPc), 1, 'last');
    if isempty(idx_poza), tr = 0; elseif idx_poza == length(t_r), tr = NaN; else, tr = t_r(idx_poza); end
    
    wyn_URM(i).model = mod_n; wyn_URM(i).metoda = met_n;
    wyn_URM(i).Kp = Kp; wyn_URM(i).Ti = Ti;
    wyn_URM(i).tr = tr; wyn_URM(i).ISE = ISE; wyn_URM(i).IAE = IAE;
    
    % Rysowanie na właściwej parze wykresów (osobno Styczna, osobno Dwupunktowy)
    if strcmp(mod_n, 'Styczna')
        figure(f_st_PV); plot(t_u, PV_u, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.2, 'DisplayName', opis);
        figure(f_st_CV); plot(t_u, CV_u, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.2, 'DisplayName', opis);
    else
        figure(f_2p_PV); plot(t_u, PV_u, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.2, 'DisplayName', opis);
        figure(f_2p_CV); plot(t_u, CV_u, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.2, 'DisplayName', opis);
    end
end

% 5. Wykończenie i eksport wszystkich 4 wykresów (300 DPI pod LaTeX)
t_sp = [0, czas_skok, czas_skok, czas];
y_sp = [Twew0, Twew0, Twew0 + dSPc, Twew0 + dSPc];

% Wykres 1: Styczna PV
figure(f_st_PV);
plot(t_sp, y_sp, 'k:', 'LineWidth', 1.8, 'DisplayName', 'Wartość zadana SP');
yline(Twew0 + 1.05*dSPc, 'k--', 'Alpha', 0.35, 'HandleVisibility', 'off');
yline(Twew0 + 0.95*dSPc, 'k--', 'Alpha', 0.35, 'HandleVisibility', 'off');
xlim([400, czas]); legend('show', 'Location', 'southeast', 'FontSize', 9.5);
exportgraphics(f_st_PV, 'URM_styczna_PV.png', 'Resolution', 300);

% Wykres 2: Styczna CV
figure(f_st_CV);
yline(Fw0, 'k:', 'Punkt pracy F_{w0}', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlim([400, czas]); legend('show', 'Location', 'northeast', 'FontSize', 9.5);
exportgraphics(f_st_CV, 'URM_styczna_CV.png', 'Resolution', 300);

% Wykres 3: Dwupunktowa PV
figure(f_2p_PV);
plot(t_sp, y_sp, 'k:', 'LineWidth', 1.8, 'DisplayName', 'Wartość zadana SP');
yline(Twew0 + 1.05*dSPc, 'k--', 'Alpha', 0.35, 'HandleVisibility', 'off');
yline(Twew0 + 0.95*dSPc, 'k--', 'Alpha', 0.35, 'HandleVisibility', 'off');
xlim([400, czas]); legend('show', 'Location', 'southeast', 'FontSize', 9.5);
exportgraphics(f_2p_PV, 'URM_dwupunktowa_PV.png', 'Resolution', 300);

% Wykres 4: Dwupunktowa CV
figure(f_2p_CV);
yline(Fw0, 'k:', 'Punkt pracy F_{w0}', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlim([400, czas]); legend('show', 'Location', 'northeast', 'FontSize', 9.5);
exportgraphics(f_2p_CV, 'URM_dwupunktowa_CV.png', 'Resolution', 300);

% 6. Tabela wyników w konsoli i gotowy kod LaTeX
if uzyto_simulinka
    fprintf('\n[OK] Dane pobrane bezposrednio z symulacji pliku model.mdl.\n');
else
    fprintf('\n[INFO] Dane wyliczone analitycznie (dodaj bloki To Workspace: PV_urm, dCV_urm w model.mdl, aby logowac z Simulinka).\n');
end

fprintf('\n=================================================================================\n');
fprintf('%-12s | %-16s | %-9s | %-8s | %-10s | %-8s | %-8s\n', ...
        'Model', 'Metoda', 'Kp', 'Ti [s]', 'tr(5%)[s]', 'ISE', 'IAE');
fprintf('---------------------------------------------------------------------------------\n');
for k = 1:n
    fprintf('%-12s | %-16s | %9.5f | %8.2f | %10.1f | %8.2f | %8.2f\n', ...
            wyn_URM(k).model, wyn_URM(k).metoda, wyn_URM(k).Kp, wyn_URM(k).Ti, ...
            wyn_URM(k).tr, wyn_URM(k).ISE, wyn_URM(k).IAE);
end
fprintf('=================================================================================\n');

fprintf('\n--- GOTOWY KOD TABELI DO LATEXA ---\n');
fprintf('\\begin{table}[H]\n\\centering\n');
fprintf('\\caption{Zestawienie nastaw regulatora PI i wskaźników jakości dla URM}\n');
fprintf('\\label{tab:urm_wyniki}\n');
fprintf('\\begin{tabular}{|l|l|c|c|c|c|c|}\n\\hline\n');
fprintf('\\textbf{Model} & \\textbf{Metoda} & \\(K_p\\) & \\(T_i\\,[\\mathrm{s}]\\) & \\(t_r\\,[\\mathrm{s}]\\) & \\textbf{ISE} & \\textbf{IAE} \\\\ \\hline\n');
for k = 1:n
    fprintf('%s & %s & %.5f & %.2f & %.1f & %.2f & %.2f \\\\ \\hline\n', ...
            wyn_URM(k).model, wyn_URM(k).metoda, wyn_URM(k).Kp, wyn_URM(k).Ti, ...
            wyn_URM(k).tr, wyn_URM(k).ISE, wyn_URM(k).IAE);
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