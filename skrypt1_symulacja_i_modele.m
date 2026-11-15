clear all; close all; clc;

% ==========================================================
% 1. PARAMETRY OBIEKTU I PUNKTU PRACY
% ==========================================================
model = 'Wezel1Pom1CvwCvg';
aparam1;

% Wartości nominalne instalacji
QgN   = 2000;         % [W]
TzewN = -20;          % [C]
TwewN = 20;           % [C]
TgzN  = 90;           % [C]
TgpN  = 70;           % [C]
TwzN  = 130;          % [C]
TwpN  = 110;          % [C]

% Parametry statyczne
FgN   = QgN / (cpw * (TgzN - TgpN));
Kcg   = QgN / (TgpN - TwewN);
Kcw   = QgN / (TwewN - TzewN);
FwN   = QgN / (TwzN - TwpN) / cpw;
Kco   = QgN / (TwpN - TgzN);
atz_c = (TwzN - TzewN) / (TwewN - TzewN);
btz_c = (TwzN - TwewN) / (TwewN - TzewN);

% Warunki początkowe
SPc0  = TwewN;
Tzew0 = TzewN;
Twz0  = atz_c * SPc0 - btz_c * Tzew0;
Fg0   = FgN;
Fw0   = FwN;
Qt0   = 0;

% Równania stanu punktu pracy
A = [(-Kcg-Kcw)/Cvws,  Kcg/Cvws,                 0,                  0;
      Kcg/Cvg,         (-cpw*Fg0-Kcg)/Cvg,       cpw*Fg0/Cvg,        0;
      0,               cpw*Fg0/Cvo,              (-cpw*Fg0-Kco)/Cvo, Kco/Cvo;
      0,               0,                        Kco/Cvo,            (-cpw*Fw0-Kco)/Cvo];
B = [0, Kcw/Cvws; 
     0, 0; 
     0, 0; 
     cpw*Fw0/Cvo, 0];
u0 = [Twz0; Tzew0];
x0 = -A^(-1) * B * u0;
Twew0 = x0(1); 
Tgp0  = x0(2); 
Tgz0  = x0(3); 
Twp0  = x0(4);

% Wymuszenie skokowe (identyfikacja w układzie otwartym)
dSPc     = 0;
dTzew    = 0;
dFg      = 0;
dQt      = 0;
dFw      = 0.5 * FwN;
czas     = 3000;
czas_skok = 500;

% Zmienne pomocnicze dla bloków w modelu Simulinka
kp = 0; Kp = 0; Ti = 1000;

% ==========================================================
% 2. SYMULACJA OBIEKTU W SIMULINKU
% ==========================================================
[t] = sim(model, czas);

y0      = aTwew(1);
y_ust   = aTwew(end);
dy_ust  = y_ust - y0;
K_stat  = dy_ust / dFw;

% ==========================================================
% A. METODA DWUPUNKTOWA (28.3% oraz 63.2%)
% ==========================================================
y_283 = y0 + 0.283 * dy_ust;
y_632 = y0 + 0.632 * dy_ust;

idx_283 = find(aTwew >= y_283, 1, 'first');
idx_632 = find(aTwew >= y_632, 1, 'first');

t_p1 = t(idx_283); y_p1 = aTwew(idx_283);
t_p2 = t(idx_632); y_p2 = aTwew(idx_632);

t1 = t_p1 - czas_skok;
t2 = t_p2 - czas_skok;

K_2p  = K_stat;
T_2p  = 1.5 * (t2 - t1);
T0_2p = t2 - T_2p;

% Odpowiedź modelu 2-punktowego
y_model_2p = zeros(size(t));
idx_akt_2p = t > (czas_skok + T0_2p);
y_model_2p(idx_akt_2p) = K_2p * dFw * (1 - exp(-(t(idx_akt_2p) - czas_skok - T0_2p) / T_2p));
Twew_model_2p = y0 + y_model_2p;

% ==========================================================
% B. METODA STYCZNEJ
% ==========================================================
dt_vec = diff(t);
dy_vec = diff(aTwew);
pochodna = dy_vec ./ dt_vec;
pochodna(dt_vec == 0) = 0;   % Zabezpieczenie przed powtórzonym czasem w chwili skoku

[max_nachylenie, idx_przeg] = max(pochodna);
t_przeg = t(idx_przeg);
y_przeg = aTwew(idx_przeg);

t_start_st  = t_przeg - (y_przeg - y0) / max_nachylenie;
t_koniec_st = t_przeg + (y_ust - y_przeg) / max_nachylenie;

K_st  = K_stat;
T0_st = t_start_st - czas_skok;
T_st  = t_koniec_st - t_start_st;

% Odpowiedź modelu ze stycznej
y_model_st = zeros(size(t));
idx_akt_st = t > (czas_skok + T0_st);
y_model_st(idx_akt_st) = K_st * dFw * (1 - exp(-(t(idx_akt_st) - czas_skok - T0_st) / T_st));
Twew_model_st = y0 + y_model_st;

% ==========================================================
% C. METODA MOMENTÓW (4 parametry: b1, b0, a2, a1)
% ==========================================================
% Selekcja próbek ściśle po chwili skoku z unikalnymi czasami
idx_m = t >= czas_skok;
t_m   = t(idx_m) - czas_skok;
dy_m  = aTwew(idx_m) - y0;
e_m   = (dy_ust - dy_m) / dy_ust;

S0 = trapz(t_m, e_m);
S1 = trapz(t_m, t_m .* e_m);
S2 = trapz(t_m, (t_m.^2 / 2) .* e_m);

b0 = K_stat;
rozw = [S0, -1; S1, -S0] \ [S1; S2];
a1 = rozw(1);
a2 = rozw(2);
b1 = (a1 - S0) * b0;

% ==========================================================
% D. TRANSMITANCJE I ZAPIS DANYCH
% ==========================================================
s = tf('s');
G_2p  = (K_2p / (T_2p * s + 1)) * exp(-T0_2p * s);
G_st  = (K_st / (T_st * s + 1)) * exp(-T0_st * s);
G_mom = (b1 * s + b0) / (a2 * s^2 + a1 * s + 1);

% Wyliczenie odpowiedzi skokowej dla metody momentów:
% 1. Generujemy wektor czasu o równomiernym kroku (wymóg funkcji step)
t_max_m = czas - czas_skok;
t_rownomierny = linspace(0, t_max_m, 3000)';

% 2. Odpowiedź skokowa na regularnej siatce
y_step_raw = dFw * step(G_mom, t_rownomierny);

% 3. Interpolacja na oryginalne punkty czasowe z Simulinka
dy_model_mom = zeros(size(t));
idx_m = t >= czas_skok;
dy_model_mom(idx_m) = interp1(t_rownomierny, y_step_raw, t(idx_m) - czas_skok, 'linear', 'extrap');
Twew_model_mom = y0 + dy_model_mom;

% Zapis wszystkich wyliczonych danych i modeli do pliku .mat
save('modele_transmitancje.mat', ...
     't', 'aTwew', 'y0', 'y_ust', 'dy_ust', 'dFw', 'czas_skok', ...
     'K_2p', 'T_2p', 'T0_2p', 't_p1', 'y_p1', 't_p2', 'y_p2', 'Twew_model_2p', 'G_2p', ...
     'K_st', 'T_st', 'T0_st', 't_przeg', 'y_przeg', 't_start_st', 't_koniec_st', 'Twew_model_st', 'G_st', ...
     'b0', 'b1', 'a1', 'a2', 'Twew_model_mom', 'G_mom');

% Wyświetlenie podsumowania w oknie Command Window
fprintf('\n================== WYNIKI IDENTYFIKACJI ==================\n');
fprintf('Wzmocnienie statyczne K = %.4f C/(kg/s)\n', K_stat);
fprintf('Metoda 2-punktowa: T = %.2f s, T0 = %.2f s\n', T_2p, T0_2p);
fprintf('Metoda stycznej:   T = %.2f s, T0 = %.2f s\n', T_st, T0_st);
fprintf('Metoda momentow:   G(s) = (%.2f*s + %.2f) / (%.2f*s^2 + %.2f*s + 1)\n', b1, b0, a2, a1);
fprintf('==========================================================\n');
fprintf('Plik "modele_transmitancje.mat" zostal pomyslnie utworzony.\n');