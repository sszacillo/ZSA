% =========================================================================
% PROJEKTOWANIE REGULATORA PI: URM i URO
% Modele: a) metoda stycznej, b) metoda dwupunktowa
% Metody: a) Ziegler-Nichols, b) CHR (20%), c) pidtune
% Wykresy: PV (Twew) i CV (Fw) + Tabela (model, nastawy, tr, ISE, IAE)
% =========================================================================
clear all; close all; clc;

% --- 1. KONFIGURACJA MODELI I PUNKTU PRACY ---
if exist('modele_transmitancje.mat', 'file')
    load('modele_transmitancje.mat');
else
    K_2p = 155.42; T_2p = 426.00; T0_2p = 96.00;
    K_st = 155.42; T_st = 618.08; T0_st = 62.15;
end

% Inicjalizacja parametrów obiektu nieliniowego (dla URO w Simulinku)
if exist('aparam1.m', 'file')
    aparam1;
    QgN = 2000; TzewN = -20; TwewN = 20;
    TgzN = 90;  TgpN = 70;   TwzN = 130; TwpN = 110;
    FgN = QgN / (cpw * (TgzN - TgpN));
    Kcg = QgN / (TgpN - TwewN);
    Kcw = QgN / (TwewN - TzewN);
    FwN = QgN / (TwzN - TwpN) / cpw;
    Kco = QgN / (TwpN - TgzN);
    atz_c = (TwzN - TzewN) / (TwewN - TzewN);
    btz_c = (TwzN - TwewN) / (TwewN - TzewN);

    SPc0 = TwewN; Tzew0 = TzewN;
    Twz0 = atz_c * SPc0 - btz_c * Tzew0;
    Fg0 = FgN; Fw0 = FwN; Qt0 = 0;

    A = [(-Kcg-Kcw)/Cvws,  Kcg/Cvws,                 0,                  0;
          Kcg/Cvg,         (-cpw*Fg0-Kcg)/Cvg,       cpw*Fg0/Cvg,        0;
          0,               cpw*Fg0/Cvo,              (-cpw*Fg0-Kco)/Cvo, Kco/Cvo;
          0,               0,                        Kco/Cvo,            (-cpw*Fw0-Kco)/Cvo];
    B = [0, Kcw/Cvws; 0, 0; 0, 0; cpw*Fw0/Cvo, 0];
    u0 = [Twz0; Tzew0];
    x0 = -A^(-1) * B * u0;
    Twew0 = x0(1); Tgp0 = x0(2); Tgz0 = x0(3); Twp0 = x0(4);
else
    Twew0 = 20.0;
    Fw0   = 0.023866;
end

% Nazwa pliku Simulinka dla URO (zmień, jeśli w repozytorium nazywa się inaczej)
model_URO = 'Wezel1Pom1CvwCvg_URO';

% Parametry eksperymentu regulacji (skok wartości zadanej SP)
czas      = 4500;  % Czas symulacji [s]
czas_skok = 500;   % Chwila skoku wartości zadanej [s]
dSPc      = 1.0;   % Skok temperatury zadanej o +1°C (20°C -> 21°C)
dTzew = 0; dFg = 0; dFw = 0; dQt = 0;
wp = Fw0;          % Warunek początkowy integratora dla URO

% --- 2. TRANSMITANCJE I OBLICZENIE NASTAW (Z-N, CHR 20%, pidtune) ---
s = tf('s');
G_st = (K_st / (T_st * s + 1)) * exp(-T0_st * s);
G_2p = (K_2p / (T_2p * s + 1)) * exp(-T0_2p * s);

% a) Na bazie modelu ZE STYCZNEJ
Kp_st_ZN  = (0.9 * T_st) / (K_st * T0_st);
Ti_st_ZN  = 3.33 * T0_st;

Kp_st_CHR = (0.6 * T_st) / (K_st * T0_st);   % CHR 20% przeregulowania (slajd 8)
Ti_st_CHR = 1.0 * T_st;
% Alternatywnie CHR 0%: Kp = (0.35*T_st)/(K_st*T0_st); Ti = 1.2*T_st;
% Alternatywnie SIMC (Tz=T0): Kp = T_st/(K_st*2*T0_st); Ti = min(T_st, 8*T0_st);

C_st_pt   = pidtune(G_st, 'PI');
Kp_st_pt  = C_st_pt.Kp;
Ti_st_pt  = C_st_pt.Kp / C_st_pt.Ki;

% b) Na bazie modelu DWUPUNKTOWEGO
Kp_2p_ZN  = (0.9 * T_2p) / (K_2p * T0_2p);
Ti_2p_ZN  = 3.33 * T0_2p;

Kp_2p_CHR = (0.6 * T_2p) / (K_2p * T0_2p);   % CHR 20% przeregulowania (slajd 8)
Ti_2p_CHR = 1.0 * T_2p;

C_2p_pt   = pidtune(G_2p, 'PI');
Kp_2p_pt  = C_2p_pt.Kp;
Ti_2p_pt  = C_2p_pt.Kp / C_2p_pt.Ki;

% Tabela konfiguracji:
% {Model, Metoda, Kp, Ti, G_model, Styl, Kolor}
konf = {
    'Styczna',     'Ziegler-Nichols', Kp_st_ZN,  Ti_st_ZN,  G_st, '--', [0.85 0.15 0.15];
    'Styczna',     'CHR (20%)',       Kp_st_CHR, Ti_st_CHR, G_st, '--', [0.00 0.45 0.74];
    'Styczna',     'pidtune',         Kp_st_pt,  Ti_st_pt,  G_st, '--', [0.15 0.65 0.20];
    'Dwupunktowy', 'Ziegler-Nichols', Kp_2p_ZN,  Ti_2p_ZN,  G_2p, '-',  [0.85 0.15 0.15];
    'Dwupunktowy', 'CHR (20%)',       Kp_2p_CHR, Ti_2p_CHR, G_2p, '-',  [0.00 0.45 0.74];
    'Dwupunktowy', 'pidtune',         Kp_2p_pt,  Ti_2p_pt,  G_2p, '-',  [0.15 0.65 0.20]
};

n = size(konf, 1);
t_urm = (0:1:czas)';
idx_skok_urm = t_urm >= czas_skok;
t_po_urm     = t_urm(idx_skok_urm) - czas_skok;
SP_urm       = Twew0 + dSPc * double(idx_skok_urm);

% Sprawdzenie dostępności modelu URO w Simulinku
ma_URO = exist([model_URO '.slx'], 'file') || exist([model_URO '.mdl'], 'file');

% Przygotowanie figur dla URM
fig_PV_urm = figure('Color', 'w', 'Position', [60, 120, 960, 480]); hold on; grid on; box on;
title('URM (Układ Regulacji Modelu) - porównanie wielkości regulowanej PV (T_{wew})');
xlabel('Czas t [s]'); ylabel('Temperatura T_{wew} [°C]');

fig_CV_urm = figure('Color', 'w', 'Position', [90, 90, 960, 480]); hold on; grid on; box on;
title('URM (Układ Regulacji Modelu) - porównanie sygnału sterującego CV (F_w)');
xlabel('Czas t [s]'); ylabel('Przepływ F_w [kg/s]');

if ma_URO
    fig_PV_uro = figure('Color', 'w', 'Position', [120, 60, 960, 480]); hold on; grid on; box on;
    title('URO (Układ Regulacji Obiektu) - porównanie wielkości regulowanej PV (T_{wew})');
    xlabel('Czas t [s]'); ylabel('Temperatura T_{wew} [°C]');

    fig_CV_uro = figure('Color', 'w', 'Position', [150, 30, 960, 480]); hold on; grid on; box on;
    title('URO (Układ Regulacji Obiektu) - porównanie sygnału sterującego CV (F_w)');
    xlabel('Czas t [s]'); ylabel('Przepływ F_w [kg/s]');
end

wyn_URM = struct();
wyn_URO = struct();

% --- 3. PĘTLA SYMULACJI DLA WSZYSTKICH 6 WARIANTÓW ---
for i = 1:n
    mod_nazwa = konf{i, 1};
    met_nazwa = konf{i, 2};
    Kp        = konf{i, 3};
    kp        = Kp;   % na wypadek małych liter w bloku Simulinka
    Ti        = konf{i, 4};
    Gm        = konf{i, 5};
    styl      = konf{i, 6};
    kolor     = konf{i, 7};
    
    etykieta = sprintf('%s | %s (K_p=%.4f, T_i=%.1f s)', mod_nazwa, met_nazwa, Kp, Ti);
    
    % --- A. SYMULACJA URM (model liniowy z opóźnieniem) ---
    R_PI   = Kp * (1 + 1 / (Ti * s));
    T_PV   = feedback(R_PI * Gm, 1);
    T_CV   = feedback(R_PI, Gm);
    
    PV_urm = Twew0 * ones(size(t_urm));
    CV_urm = Fw0   * ones(size(t_urm));
    PV_urm(idx_skok_urm) = Twew0 + dSPc * step(T_PV, t_po_urm);
    CV_urm(idx_skok_urm) = Fw0   + dSPc * step(T_CV, t_po_urm);
    
    [tr_u, ise_u, iae_u] = oblicz_wskazniki(t_urm, PV_urm, czas_skok, Twew0 + dSPc, dSPc, 0.05);
    wyn_URM(i).model = mod_nazwa; wyn_URM(i).metoda = met_nazwa;
    wyn_URM(i).Kp = Kp; wyn_URM(i).Ti = Ti;
    wyn_URM(i).tr = tr_u; wyn_URM(i).ISE = ise_u; wyn_URM(i).IAE = iae_u;
    
    figure(fig_PV_urm);
    plot(t_urm, PV_urm, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.0, 'DisplayName', etykieta);
    figure(fig_CV_urm);
    plot(t_urm, CV_urm, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.0, 'DisplayName', etykieta);
    
    % --- B. SYMULACJA URO (obiekt nieliniowy w Simulinku) ---
    if ma_URO
        [t_sim] = sim(model_URO, czas);
        PV_uro = aTwew;
        CV_uro = aFw;
        
        [tr_o, ise_o, iae_o] = oblicz_wskazniki(t_sim, PV_uro, czas_skok, Twew0 + dSPc, dSPc, 0.05);
        wyn_URO(i).model = mod_nazwa; wyn_URO(i).metoda = met_nazwa;
        wyn_URO(i).Kp = Kp; wyn_URO(i).Ti = Ti;
        wyn_URO(i).tr = tr_o; wyn_URO(i).ISE = ise_o; wyn_URO(i).IAE = iae_o;
        
        figure(fig_PV_uro);
        plot(t_sim, PV_uro, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.0, 'DisplayName', etykieta);
        figure(fig_CV_uro);
        plot(t_sim, CV_uro, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.0, 'DisplayName', etykieta);
    end
end

% --- 4. FORMATOWANIE I EKSPORT WYKRESÓW (300 DPI DO LATEXA) ---
figure(fig_PV_urm);
plot(t_urm, SP_urm, 'k:', 'LineWidth', 1.8, 'DisplayName', 'Wartość zadana SP');
yline(Twew0 + 1.05*dSPc, 'k--', 'Alpha', 0.35, 'HandleVisibility', 'off');
yline(Twew0 + 0.95*dSPc, 'k--', 'Alpha', 0.35, 'HandleVisibility', 'off');
xlim([400, czas]); legend('show', 'Location', 'southeast', 'FontSize', 9);
exportgraphics(fig_PV_urm, 'URM_porownanie_PV.png', 'Resolution', 300);

figure(fig_CV_urm);
yline(Fw0, 'k:', 'HandleVisibility', 'off');
xlim([400, czas]); legend('show', 'Location', 'northeast', 'FontSize', 9);
exportgraphics(fig_CV_urm, 'URM_porownanie_CV.png', 'Resolution', 300);

if ma_URO
    figure(fig_PV_uro);
    plot(t_urm, SP_urm, 'k:', 'LineWidth', 1.8, 'DisplayName', 'Wartość zadana SP');
    yline(Twew0 + 1.05*dSPc, 'k--', 'Alpha', 0.35, 'HandleVisibility', 'off');
    yline(Twew0 + 0.95*dSPc, 'k--', 'Alpha', 0.35, 'HandleVisibility', 'off');
    xlim([400, czas]); legend('show', 'Location', 'southeast', 'FontSize', 9);
    exportgraphics(fig_PV_uro, 'URO_porownanie_PV.png', 'Resolution', 300);

    figure(fig_CV_uro);
    yline(Fw0, 'k:', 'HandleVisibility', 'off');
    xlim([400, czas]); legend('show', 'Location', 'northeast', 'FontSize', 9);
    exportgraphics(fig_CV_uro, 'URO_porownanie_CV.png', 'Resolution', 300);
else
    fprintf('\n[INFO] Nie znaleziono modelu "%s.slx" - wygenerowano wyniki dla URM.\n', model_URO);
    fprintf('       Zapisz swoj schemat URO pod ta nazwa i uruchom skrypt ponownie.\n');
end

% --- 5. WYPISANIE TABEL W KONSOLI I KODU LATEX ---
wypisz_tabele('URM (Układ Regulacji na Modelu)', wyn_URM);
if ma_URO
    wypisz_tabele('URO (Układ Regulacji na Obiekcie)', wyn_URO);
end

% =========================================================================
% FUNKCJE POMOCNICZE
% =========================================================================
function [tr, ISE, IAE] = oblicz_wskazniki(t, y, t_skok, SP_doc, dSP, eps_proc)
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

function wypisz_tabele(tytul, W)
    fprintf('\n==================================================================================\n');
    fprintf('  TABELA WYNIKÓW: %s\n', tytul);
    fprintf('==================================================================================\n');
    fprintf('%-12s | %-16s | %-9s | %-8s | %-10s | %-8s | %-8s\n', ...
            'Model', 'Metoda', 'Kp', 'Ti [s]', 'tr(5%)[s]', 'ISE', 'IAE');
    fprintf('----------------------------------------------------------------------------------\n');
    for k = 1:length(W)
        fprintf('%-12s | %-16s | %9.5f | %8.2f | %10.1f | %8.2f | %8.2f\n', ...
                W(k).model, W(k).metoda, W(k).Kp, W(k).Ti, W(k).tr, W(k).ISE, W(k).IAE);
    end
    
    fprintf('\n--- KOD LATEX (%s) ---\n', tytul);
    fprintf('\\begin{table}[H]\n\\centering\n');
    fprintf('\\caption{Zestawienie nastaw i wskaźników jakości regulacji -- %s}\n', tytul);
    fprintf('\\begin{tabular}{|l|l|c|c|c|c|c|}\n\\hline\n');
    fprintf('\\textbf{Model} & \\textbf{Metoda} & \\(K_p\\) & \\(T_i\\,[\\mathrm{s}]\\) & \\(t_r\\,[\\mathrm{s}]\\) & \\textbf{ISE} & \\textbf{IAE} \\\\ \\hline\n');
    for k = 1:length(W)
        fprintf('%s & %s & %.5f & %.2f & %.1f & %.2f & %.2f \\\\ \\hline\n', ...
                W(k).model, W(k).metoda, W(k).Kp, W(k).Ti, W(k).tr, W(k).ISE, W(k).IAE);
    end
    fprintf('\\end{tabular}\n\\end{table}\n');
end