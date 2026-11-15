% =========================================================================
% KOMPLETNY SKRYPT: REGULACJA PI (URM i URO), WYKRESY PV/CV ORAZ TABELA ISE
% =========================================================================
clear all; close all; clc;

% 1. WCZYTANIE MODELI Z IDENTYFIKACJI I INICJALIZACJA PUNKTU PRACY OBIEKTU
if exist('modele_transmitancje.mat', 'file')
    load('modele_transmitancje.mat');
else
    % Wartości wyznaczone z Twojej identyfikacji (jeśli plik .mat nie jest w folderze)
    K_2p = 155.42; T_2p = 426.00; T0_2p = 96.00;
    K_st = 155.42; T_st = 618.08; T0_st = 62.15;
end

% Inicjalizacja parametrów fizycznych obiektu (do symulacji URO w Simulinku)
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

SPc0  = TwewN;
Tzew0 = TzewN;
Twz0  = atz_c * SPc0 - btz_c * Tzew0;
Fg0   = FgN;
Fw0   = FwN;
Qt0   = 0;

A = [(-Kcg-Kcw)/Cvws,  Kcg/Cvws,                 0,                  0;
      Kcg/Cvg,         (-cpw*Fg0-Kcg)/Cvg,       cpw*Fg0/Cvg,        0;
      0,               cpw*Fg0/Cvo,              (-cpw*Fg0-Kco)/Cvo, Kco/Cvo;
      0,               0,                        Kco/Cvo,            (-cpw*Fw0-Kco)/Cvo];
B = [0, Kcw/Cvws; 0, 0; 0, 0; cpw*Fw0/Cvo, 0];
u0 = [Twz0; Tzew0];
x0 = -A^(-1) * B * u0;
Twew0 = x0(1); Tgp0 = x0(2); Tgz0 = x0(3); Twp0 = x0(4);

% Parametry skoku wartości zadanej SP (wymuszenie w układzie zamkniętym)
czas      = 4500;   % Czas symulacji [s]
czas_skok = 500;    % Chwila skoku wartości zadanej [s]
dSPc      = 1.0;    % Skok wartości zadanej o +1°C (z 20°C na 21°C)
dTzew     = 0; dFg = 0; dFw = 0; dQt = 0;

% =========================================================================
% 2. DOBÓR NASTAW REGULATORA PI (Ziegler-Nichols, CHR, pidtune)
% =========================================================================
s = tf('s');
G_2p_tf = (K_2p / (T_2p * s + 1)) * exp(-T0_2p * s);
G_st_tf = (K_st / (T_st * s + 1)) * exp(-T0_st * s);

% --- A. Na bazie modelu DWUPUNKTOWEGO (K=155.42, T=426.00, T0=96.00) ---
% 1) Metoda Zieglera-Nicholsa (Z-N)
Kp_ZN_2p  = (0.9 * T_2p) / (K_2p * T0_2p);
Ti_ZN_2p  = 3.33 * T0_2p;

% 2) Metoda CHR (Chien-Hrones-Reswick, dla 20% przeregulowania wartości zadanej)
Kp_CHR_2p = (0.6 * T_2p) / (K_2p * T0_2p);
Ti_CHR_2p = 1.0 * T_2p;

% 3) Opcja pidtune
[C_pt_2p, ~] = pidtune(G_2p_tf, 'PI');
Kp_pt_2p  = C_pt_2p.Kp;
Ti_pt_2p  = C_pt_2p.Kp / C_pt_2p.Ki;

% --- B. Na bazie modelu ZE STYCZNEJ (K=155.42, T=618.08, T0=62.15) ---
% 1) Metoda Zieglera-Nicholsa (Z-N)
Kp_ZN_st  = (0.9 * T_st) / (K_st * T0_st);
Ti_ZN_st  = 3.33 * T0_st;

% 2) Metoda CHR (20% przeregulowania)
Kp_CHR_st = (0.6 * T_st) / (K_st * T0_st);
Ti_CHR_st = 1.0 * T_st;

% 3) Opcja pidtune
[C_pt_st, ~] = pidtune(G_st_tf, 'PI');
Kp_pt_st  = C_pt_st.Kp;
Ti_pt_st  = C_pt_st.Kp / C_pt_st.Ki;

% Zestawienie wszystkich przypadków w tablicy komórkowej:
% {Etykieta, Model_bazowy, Metoda, Kp, Ti, Obiekt_G_URM, Styl_linii, Kolor}
przypadki = {
    '2-pkt: Z-N',      'Dwupunktowy', 'Ziegler-Nichols', Kp_ZN_2p,  Ti_ZN_2p,  G_2p_tf, '-',  [0.85 0.10 0.10];
    '2-pkt: CHR',      'Dwupunktowy', 'CHR (20%)',       Kp_CHR_2p, Ti_CHR_2p, G_2p_tf, '-',  [0.00 0.45 0.74];
    '2-pkt: pidtune',  'Dwupunktowy', 'pidtune',         Kp_pt_2p,  Ti_pt_2p,  G_2p_tf, '-',  [0.10 0.65 0.15];
    'Styczna: Z-N',    'Styczna',     'Ziegler-Nichols', Kp_ZN_st,  Ti_ZN_st,  G_st_tf, '--', [0.90 0.40 0.00];
    'Styczna: CHR',    'Styczna',     'CHR (20%)',       Kp_CHR_st, Ti_CHR_st, G_st_tf, '--', [0.49 0.18 0.56];
    'Styczna: pidtune','Styczna',     'pidtune',         Kp_pt_st,  Ti_pt_st,  G_st_tf, '--', [0.00 0.65 0.65]
};

% =========================================================================
% 3. SYMULACJE URM (ANALITYCZNE) ORAZ URO (SIMULINK)
% =========================================================================
% Sprawdzenie, pod jaką nazwą zapisany jest model URO w Simulinku
if exist('Wezel1Pom1CvwCvg_URO.slx', 'file')
    model_URO = 'Wezel1Pom1CvwCvg_URO';
elseif exist('Wezel1Pom1CvwCvg.slx', 'file')
    model_URO = 'Wezel1Pom1CvwCvg';
else
    model_URO = '';
end

% Przygotowanie wspólnego wektora czasu do symulacji URM
t_vec = linspace(0, czas, 4501)';
idx_po_skoku = t_vec >= czas_skok;
t_po_skoku   = t_vec(idx_po_skoku) - czas_skok;
SP_vec       = Twew0 + dSPc * double(t_vec >= czas_skok);

% Przygotowanie okien wykresów
% Wykres 1: Porównanie PV i CV dla wszystkich nastaw (URO / URM)
fig_PV = figure('Color', 'w', 'Position', [80, 80, 980, 500]);
hold on; grid on; box on;
title('Porównanie wartości regulowanej PV (Temperatura wewnętrzna T_{wew})', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11);
ylabel('Temperatura T_{wew} [°C]', 'FontSize', 11);

fig_CV = figure('Color', 'w', 'Position', [120, 120, 980, 500]);
hold on; grid on; box on;
title('Porównanie sygnału sterującego CV (Przepływ wody sieciowej F_w)', 'FontSize', 12);
xlabel('Czas t [s]', 'FontSize', 11);
ylabel('Przepływ F_w [kg/s]', 'FontSize', 11);

% Tablica na wyniki do tabeli końcowej
liczba_przyp = size(przypadki, 1);
wyniki_tr  = zeros(liczba_przyp, 1);
wyniki_ISE = zeros(liczba_przyp, 1);
wyniki_IAE = zeros(liczba_przyp, 1);

% Zapis przebiegów URM vs URO dla wybranego reprezentatywnego przypadku (np. 2-pkt Z-N i CHR)
PV_URM_zapis = zeros(length(t_vec), liczba_przyp);
CV_URM_zapis = zeros(length(t_vec), liczba_przyp);
PV_URO_zapis = cell(liczba_przyp, 1);
CV_URO_zapis = cell(liczba_przyp, 1);
t_URO_zapis  = cell(liczba_przyp, 1);

czy_URO_dziala = false;

for i = 1:liczba_przyp
    etykieta = przypadki{i, 1};
    Kp = przypadki{i, 4};
    kp = Kp; % dla kompatybilności z małymi i wielkimi literami w Simulinku
    Ti = przypadki{i, 5};
    G_aktualne = przypadki{i, 6};
    styl  = przypadki{i, 7};
    kolor = przypadki{i, 8};
    
    % --- A. Symulacja URM (Układ Regulacji Modelu) ---
    % Transmitancja regulatora PI-ISA: R(s) = Kp * (1 + 1/(Ti*s))
    R_PI = Kp * (1 + 1 / (Ti * s));
    
    % Transmitancja zamkniętego układu URM dla PV (od dSP do dPV) i CV (od dSP do dCV)
    T_URM_PV = feedback(R_PI * G_aktualne, 1);
    T_URM_CV = feedback(R_PI, G_aktualne);
    
    dPV_urm = dSPc * step(T_URM_PV, t_po_skoku);
    dCV_urm = dSPc * step(T_URM_CV, t_po_skoku);
    
    PV_urm_pelne = Twew0 * ones(size(t_vec));
    CV_urm_pelne = Fw0   * ones(size(t_vec));
    PV_urm_pelne(idx_po_skoku) = Twew0 + dPV_urm;
    CV_urm_pelne(idx_po_skoku) = Fw0   + dCV_urm;
    
    PV_URM_zapis(:, i) = PV_urm_pelne;
    CV_URM_zapis(:, i) = CV_urm_pelne;
    
    % --- B. Symulacja URO w Simulinku (jeśli model ma podpięty regulator) ---
    if ~isempty(model_URO)
        try
            [t_sim] = sim(model_URO, czas);
            % Sprawdzenie, czy układ zareagował na skok dSPc (czy regulator jest wpięty)
            if max(abs(aTwew - Twew0)) > 0.05
                czy_URO_dziala = true;
                t_biez  = t_sim;
                PV_biez = aTwew;
                CV_biez = aFw;
                t_URO_zapis{i}  = t_biez;
                PV_URO_zapis{i} = PV_biez;
                CV_URO_zapis{i} = CV_biez;
            end
        catch
            czy_URO_dziala = false;
        end
    end
    
    % Wybór sygnału do głównej analizy wskaźników jakości (URO jeśli dostępne, inaczej URM)
    if czy_URO_dziala
        t_analiza  = t_URO_zapis{i};
        PV_analiza = PV_URO_zapis{i};
        CV_analiza = CV_URO_zapis{i};
    else
        t_analiza  = t_vec;
        PV_analiza = PV_urm_pelne;
        CV_analiza = CV_urm_pelne;
    end
    
    % --- C. Obliczenie wskaźników jakości regulacji (tr, ISE, IAE) ---
    idx_reg = t_analiza >= czas_skok;
    t_reg   = t_analiza(idx_reg) - czas_skok;
    SP_reg  = Twew0 + dSPc;
    e_reg   = SP_reg - PV_analiza(idx_reg);
    
    % Całkowe wskaźniki jakości (metoda trapezów)
    ISE = trapz(t_reg, e_reg.^2);
    IAE = trapz(t_reg, abs(e_reg));
    
    % Czas regulacji tr (wejście w pasmo +/- 5% wartości skoku dSPc)
    pasmo_5proc = 0.05 * abs(dSPc);
    idx_poza = find(abs(e_reg) > pasmo_5proc, 1, 'last');
    if isempty(idx_poza)
        tr = 0;
    elseif idx_poza == length(t_reg)
        tr = NaN; % nie ustalił się w czasie symulacji
    else
        tr = t_reg(idx_poza);
    end
    
    wyniki_tr(i)  = tr;
    wyniki_ISE(i) = ISE;
    wyniki_IAE(i) = IAE;
    
    % --- D. Rysowanie przebiegów na głównych wykresach PV i CV ---
    opis_leg = sprintf('%s (K_p=%.4f, T_i=%.1f s)', etykieta, Kp, Ti);
    
    figure(fig_PV);
    plot(t_analiza, PV_analiza, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.0, 'DisplayName', opis_leg);
    
    figure(fig_CV);
    plot(t_analiza, CV_analiza, 'LineStyle', styl, 'Color', kolor, 'LineWidth', 2.0, 'DisplayName', opis_leg);
end

% Dokończenie wykresu PV (dodanie linii wartości zadanej SP i pasma +/- 5%)
figure(fig_PV);
plot(t_vec, SP_vec, 'k:', 'LineWidth', 1.8, 'DisplayName', 'Wartość zadana SP');
yline(Twew0 + dSPc * 1.05, 'k--', 'Pasmo +5%', 'Alpha', 0.4, 'HandleVisibility', 'off');
yline(Twew0 + dSPc * 0.95, 'k--', 'Pasmo -5%', 'Alpha', 0.4, 'HandleVisibility', 'off');
xlim([400, czas]);
ylim([Twew0 - 0.1, max(PV_URM_zapis(:)) + 0.2]);
legend('show', 'Location', 'southeast', 'FontSize', 9.5);
exportgraphics(fig_PV, 'regulacja_porownanie_PV.png', 'Resolution', 300);

% Dokończenie wykresu CV
figure(fig_CV);
yline(Fw0, 'k:', 'Punkt pracy F_{w0}', 'HandleVisibility', 'off');
xlim([400, czas]);
legend('show', 'Location', 'northeast', 'FontSize', 9.5);
exportgraphics(fig_CV, 'regulacja_porownanie_CV.png', 'Resolution', 300);

% =========================================================================
% 4. DODATKOWY WYKRES: PORÓWNANIE URM (Model) vs URO (Obiekt) - jeśli URO aktywne
% =========================================================================
if czy_URO_dziala
    fig_urm_uro = figure('Color', 'w', 'Position', [150, 150, 980, 500]);
    hold on; grid on; box on;
    title('Porównanie odpowiedzi URM (transmitancja) i URO (obiekt nieliniowy) dla metody Z-N', 'FontSize', 12);
    xlabel('Czas t [s]', 'FontSize', 11);
    ylabel('Temperatura T_{wew} [°C]', 'FontSize', 11);
    
    plot(t_vec, SP_vec, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Wartość zadana SP');
    plot(t_vec, PV_URM_zapis(:, 1), 'r--', 'LineWidth', 2.0, 'DisplayName', 'URM (Model 2-pkt, Z-N)');
    plot(t_URO_zapis{1}, PV_URO_zapis{1}, 'r-', 'LineWidth', 2.2, 'DisplayName', 'URO (Obiekt, 2-pkt Z-N)');
    plot(t_vec, PV_URM_zapis(:, 4), 'b--', 'LineWidth', 2.0, 'DisplayName', 'URM (Model ze stycznej, Z-N)');
    plot(t_URO_zapis{4}, PV_URO_zapis{4}, 'b-', 'LineWidth', 2.2, 'DisplayName', 'URO (Obiekt, styczna Z-N)');
    
    xlim([400, czas]);
    legend('show', 'Location', 'southeast', 'FontSize', 10);
    exportgraphics(fig_urm_uro, 'porownanie_URM_vs_URO.png', 'Resolution', 300);
end

% =========================================================================
% 5. WYPISANIE TABELI WYNIKÓW W KONSOLI ORAZ KODU TABELI DO LATEXA
% =========================================================================
fprintf('\n=========================================================================================\n');
if czy_URO_dziala
    fprintf('            TABELA WYNIKÓW REGULACJI (Symulacja na obiekcie URO)\n');
else
    fprintf('            TABELA WYNIKÓW REGULACJI (Symulacja na modelu URM)\n');
end
fprintf('=========================================================================================\n');
fprintf('%-13s | %-16s | %-12s | %-9s | %-11s | %-8s | %-8s\n', ...
        'Model bazowy', 'Metoda doboru', 'Kp [kg/s*C]', 'Ti [s]', 'tr (5%) [s]', 'ISE', 'IAE');
fprintf('-----------------------------------------------------------------------------------------\n');

for i = 1:liczba_przyp
    fprintf('%-13s | %-16s | %12.5f | %9.2f | %11.1f | %8.2f | %8.2f\n', ...
            przypadki{i, 2}, przypadki{i, 3}, przypadki{i, 4}, przypadki{i, 5}, ...
            wyniki_tr(i), wyniki_ISE(i), wyniki_IAE(i));
end
fprintf('=========================================================================================\n');

% Wygenerowanie gotowego kodu tabeli LaTeX do skopiowania prosto do sprawozdania
fprintf('\n--- GOTOWY KOD TABELI DO WKLEJENIA W LATEXU ---\n\n');
fprintf('\\begin{table}[H]\n\\centering\n');
fprintf('\\caption{Zestawienie nastaw regulatora PI oraz wskaźników jakości regulacji}\n');
fprintf('\\label{tab:wyniki_regulacji}\n');
fprintf('\\begin{tabular}{|l|l|c|c|c|c|c|}\n\\hline\n');
fprintf('\\textbf{Model bazowy} & \\textbf{Metoda doboru} & \\(K_p\\,\\left[\\frac{\\mathrm{kg}}{\\mathrm{s\\cdot{}^\\circ C}}\\right]\\) & \\(T_i\\,[\\mathrm{s}]\\) & \\(t_r\\,[\\mathrm{s}]\\) & \\textbf{ISE} & \\textbf{IAE} \\\\ \\hline\n');

for i = 1:liczba_przyp
    fprintf('%s & %s & %.5f & %.2f & %.1f & %.2f & %.2f \\\\ \\hline\n', ...
            strrep(sprintf('%.5f', przypadki{i,4}), '.', '{,}'), ... % pomocnicze formatowanie
            przypadki{i, 2}, przypadki{i, 3}, przypadki{i, 4}, przypadki{i, 5}, ...
            wyniki_tr(i), wyniki_ISE(i), wyniki_IAE(i));
end
fprintf('\\end{tabular}\n\\end{table}\n');