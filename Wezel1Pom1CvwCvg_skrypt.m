
clear all; close all;
%demo dla pomieszczenia z grzejnikem wodnym
model='Wezel1Pom1CvwCvg';
aparam1;
%wartoœci nominalne
%parametry instalacji w warunkach obliczeniowych
QgN = 2000;         %[W] (zapotrezbowanie budynku)
TzewN = -20;		%[C] (temperatura na zewn¹trz)
TwewN = 20;			%[C] (zadana temperatura wewnetrzna)
TgzN = 90;              %[C]  
TgpN = 70;              %[C] 
TwzN = 130;		%[C]  
TwpN = 110;		%[C] 

%parametry statyczne
FgN = QgN/(cpw*(TgzN-TgpN));	%[kg/s]
Kcg = QgN/(TgpN-TwewN);			  %[W/K]
Kcw = QgN/(TwewN-TzewN);	    %[W/K]
FwN = QgN / (TwzN - TwpN) / cpw;
Kco = QgN / (TwpN - TgzN);          %[W/K]
atz_c = (TwzN-TzewN)/(TwewN-TzewN); btz_c = (TwzN-TwewN)/(TwewN-TzewN);
%parametry "dynamiczne" - aparam1
%==========================
%warunki pocz¹tkowe
SPc0 = TwewN;
Tzew0 = TzewN+0;	%+1
Twz0 = atz_c*SPc0 - btz_c*Tzew0;
Fg0 = FgN;
Fw0 = FwN;
Qt0 = 0;
%stan równowagi
a = cpw*Fg0; b = cpw*Fw0;

%do równañ stanu i r.statycznych X=Tgp, Twew
%u=[Twz; Tzew];	x=[Twew; Tgp; Tgz; Twp], y=x
A = [(-Kcg-Kcw)/Cvws, Kcg/Cvws,                 0,                                0;...
      Kcg/Cvg,                (-cpw*Fg0-Kcg)/Cvg, cpw*Fg0/Cvg,            0;...
	  0,                           cpw*Fg0/Cvo,             (-cpw*Fg0-Kco)/Cvo, Kco/Cvo; ...
  	  0,                           0,                                Kco/Cvo,                     (-cpw*Fw0-Kco)/Cvo; ];
B = [0, Kcw/Cvws; 
       0, 0;  
       0, 0;         
       cpw*Fw0/Cvo,        0];
C = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]; D=[0,0;0,0;0,0;0,0];

u0 = [Twz0; Tzew0];					%u=[Tgz; Tz],
x0 = -A^(-1)*B*u0;
Twew0 = x0(1);Tgp0 = x0(2);Tgz0 = x0(3);Twp0 = x0(4);			%x=[Tw; Tgp; Tgz; Twp], 
%==========================
%zak³ócenie
dSPc = 0;
dTzew = 0;
dFg = 0*FgN;
dFw = 0.5*FwN;

dQt = 0;

kp = (0.9 * 426) / (155.4 * 96);   % Wynik: ok. 0.0257
Ti = 3.33 * 96;
%==========================
%symulacja
czas =3000;	czas_skok = 500;
%opcje = simget(model); %opcje = simset('MaxStep', tmax, 'RelTol',terr);
[t] = sim(model, czas); 
figure(1); hold on; grid on; title('Twew'); plot(t,aTwew,'m'),
figure(2); hold on; grid on; title('Twz, Twp'); plot(t,aTwz,'r'), plot(t,aTwp,'b'), legend('Twz', 'Twp');
figure(3); hold on; grid on; title('Twew od Fw '); plot(aFw,aTwew,'m'),

%rysowanie wykresu
figure(4); 
hold on; 
grid on; 
title('OdpowiedŸ skokowa Twew na skok przep³ywu dFw = 0.5*FwN'); 

% Rysowanie wykresów: skok przerywan¹ czerwon¹ lini¹ ('r--'), odpowiedŸ ci¹g³¹ niebiesk¹ ('b')
h1 = plot(t, TwewN+aFw, 'r--', 'LineWidth', 2);
h2 = plot(t, aTwew, 'b', 'LineWidth', 2);

% Podpisy osi
xlabel('Czas t [s]');
ylabel('Temperatura Twew [°C] / Przep³yw Fw + offset');

% ==========================================================
% IDENTYFIKACJA: METODA DWUPUNKTOWA (punkty 28.3% oraz 63.2%)
% ==========================================================
% 1. Przesuniêcie danych wzglêdem stanu pocz¹tkowego (przed skokiem)
y0 = aTwew(1);                 % Temperatura pocz¹tkowa (20 C)
y_ust = aTwew(end);            % Temperatura ustalona po skoku
dy_ust = y_ust - y0;           % Ca³kowity przyrost temperatury

% 2. Wzmocnienie statyczne K
K_2p = dy_ust / dFw;           % dFw = 0.5 * FwN

% 3. Szukanie czasów dla 28.3% oraz 63.2% przyrostu ustalonego
y_283 = y0 + 0.283 * dy_ust;
y_632 = y0 + 0.632 * dy_ust;

% Znalezienie indeksów, w których temperatura przekracza te progi
idx_283 = find(aTwew >= y_283, 1, 'first');
idx_632 = find(aTwew >= y_632, 1, 'first');

t_p1 = t(idx_283);             % rzeczywisty czas na wykresie dla 28.3%
t_p2 = t(idx_632);             % rzeczywisty czas na wykresie dla 63.2%
y_p1 = aTwew(idx_283);         % rzeczywista wartoœæ Twew dla 28.3%
y_p2 = aTwew(idx_632);         % rzeczywista wartoœæ Twew dla 63.2%

% Czas mierzony od momentu zadania skoku (czas_skok = 500 s)
t1 = t_p1 - czas_skok;         % czas osi¹gniêcia 28.3% od skoku
t2 = t_p2 - czas_skok;         % czas osi¹gniêcia 63.2% od skoku

% 4. Obliczenie sta³ej czasowej T i opóŸnienia T0
T_2p = 1.5 * (t2 - t1);
T0_2p = t2 - T_2p;

% Wyœwietlenie wyników w Command Window (z wiêksz¹ dok³adnoœci¹ dla dFw)
fprintf('\n--- WYNIKI: METODA DWUPUNKTOWA ---\n');
fprintf('Przyrost wyjœcia dY = %.4f C (dla dFw = %.4f kg/s)\n', dy_ust, dFw);
fprintf('Wzmocnienie K  = %.4f C/(kg/s)\n', K_2p);
fprintf('Czas t1 (28.3%%) = %.2f s, Czas t2 (63.2%%) = %.2f s\n', t1, t2);
fprintf('Sta³a czasowa T = %.2f s\n', T_2p);
fprintf('OpóŸnienie T0   = %.2f s\n', T0_2p);
fprintf('Transmitancja G(s) = %.2f / (%.2f*s + 1) * exp(-%.2f*s)\n', K_2p, T_2p, T0_2p);

% 5. Dorysowanie punktów pomiarowych, linii pomocniczych i wartoœci na figure(4)
figure(4);

% Linie pomocnicze rzutuj¹ce punkty na osie (szare kropkowane)
plot([t_p1, t_p1], [y0, y_p1], 'k:', 'LineWidth', 1);
plot([0, t_p1], [y_p1, y_p1], 'k:', 'LineWidth', 1);
plot([t_p2, t_p2], [y0, y_p2], 'k:', 'LineWidth', 1);
plot([0, t_p2], [y_p2, y_p2], 'k:', 'LineWidth', 1);

% Punkty 28.3% oraz 63.2%
h3 = plot(t_p1, y_p1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
h4 = plot(t_p2, y_p2, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');

% Etykiety tekstowe z wartoœciami obok punktów na wykresie
opis_p1 = sprintf('  28.3%%: (t = %.1f s, T_{wew} = %.2f °C)', t_p1, y_p1);
opis_p2 = sprintf('  63.2%%: (t = %.1f s, T_{wew} = %.2f °C)', t_p2, y_p2);

text(t_p1 + 30, y_p1 - 0.05, opis_p1, 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0 0.5 0]);
text(t_p2 + 30, y_p2 - 0.05, opis_p2, 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.6 0.5 0]);

% Legenda przypisana tylko do g³ównych krzywych i punktów (pomija linie pomocnicze)
legend([h1, h2, h3, h4], {'aFw (skok 0.5*FwN)', 'Twew', 'Punkt 28.3%', 'Punkt 63.2%'}, ...
       'FontSize', 11, 'Location', 'southeast');


% ==========================================================
% WERYFIKACJA MODELU (Porównanie obiektu z transmitancj¹)
% ==========================================================
% 1. Obliczenie odpowiedzi skokowej modelu FOPDT dla czasu t
% Wzór na odpowiedŸ skokow¹ inercji I rzêdu z opóŸnieniem:
% dy(t) = K * dFw * (1 - exp(-(t - czas_skok - T0)/T)) dla t > czas_skok + T0

y_model_2p = zeros(size(t)); % wektor na przyrosty modelu

% Indeksy czasu po uwzglêdnieniu chwili skoku i opóŸnienia T0
idx_aktywne = t > (czas_skok + T0_2p);
t_efektywny = t(idx_aktywne) - czas_skok - T0_2p;

% Obliczenie przyrostu temperatury modelu
y_model_2p(idx_aktywne) = K_2p * dFw * (1 - exp(-t_efektywny / T_2p));

% Dodanie punktu pracy (temperatury pocz¹tkowej y0 = 20 C)
Twew_model_2p = y0 + y_model_2p;

% 2. Dorysowanie odpowiedzi modelu na figure(4)
figure(4);
h5 = plot(t, Twew_model_2p, 'm-.', 'LineWidth', 2);

% Aktualizacja legendy
legend([h1, h2, h3, h4, h5], ...
       {'aFw (skok 0.5*FwN)', 'Twew (Obiekt)', 'Punkt 28.3%', 'Punkt 63.2%', 'Model 2-punktowy'}, ...
       'FontSize', 11, 'Location', 'southeast');


% ==========================================================
% IDENTYFIKACJA: METODA STYCZNEJ (w punkcie przegiêcia)
% ==========================================================
% 1. Obliczenie pochodnej dTwew/dt i znalezienie maksimum (punktu przegiêcia)
dt_vec = diff(t);
dy_vec = diff(aTwew);
pochodna = dy_vec ./ dt_vec;

[max_nachylenie, idx_przeg] = max(pochodna);
t_przeg = t(idx_przeg);
y_przeg = aTwew(idx_przeg);

% 2. Wyznaczenie przeciêcia stycznej z poziomem pocz¹tkowym y0 i ustalonym y_ust
% Równanie stycznej: y(t) = max_nachylenie * (t - t_przeg) + y_przeg
t_start_st = t_przeg - (y_przeg - y0) / max_nachylenie;       % przeciêcie z 20 C
t_koniec_st = t_przeg + (y_ust - y_przeg) / max_nachylenie;   % przeciêcie z y_ust

T0_st = t_start_st - czas_skok;
T_st = t_koniec_st - t_start_st;
K_st = K_2p; % wzmocnienie statyczne jest takie samo

fprintf('\n--- WYNIKI: METODA STYCZNEJ ---\n');
fprintf('Punkt przegiêcia: t = %.2f s, Twew = %.2f C\n', t_przeg, y_przeg);
fprintf('Sta³a czasowa T = %.2f s, OpóŸnienie T0 = %.2f s\n', T_st, T0_st);
fprintf('Transmitancja G_st(s) = %.2f / (%.2f*s + 1) * exp(-%.2f*s)\n', K_st, T_st, T0_st);

% 3. OdpowiedŸ skokowa modelu z metody stycznej (weryfikacja)
y_model_st = zeros(size(t));
idx_akt_st = t > (czas_skok + T0_st);
y_model_st(idx_akt_st) = K_st * dFw * (1 - exp(-(t(idx_akt_st) - czas_skok - T0_st) / T_st));
Twew_model_st = y0 + y_model_st;

% Dorysowanie modelu ze stycznej na wykresie figure(4)
figure(4);
h6 = plot(t, Twew_model_st, 'g--', 'LineWidth', 2);
legend([h1, h2, h3, h4, h5, h6], ...
       {'aFw (skok 0.5*FwN)', 'Twew (Obiekt)', 'Punkt 28.3%', 'Punkt 63.2%', ...
        'Model 2-punktowy', 'Model z metody stycznej'}, ...
       'FontSize', 10, 'Location', 'southeast');



Kp_st = (0.9 * T_st) / (K_st * T0_st);
Ti_st = 3.33 * T0_st;

% Przypisanie do zmiennych globalnych, jeœli u¿ywasz ich w Simulinku


Kp = (0.9 * 426) / (155.4 * 96); % Wynik: ok. 0.0257
Ti = 3.33 * 96; 

% =========================================================================
% GRAFICZNE PRZEDSTAWIENIE METODY STYCZNEJ (Osobny wykres do sprawozdania)
% =========================================================================
figure(5); hold on; grid on;
title('Identyfikacja metod¹ stycznej - konstrukcja geometryczna');
xlabel('Czas t [s]'); ylabel('Temperatura T_{wew} [°C]');

% 1. OdpowiedŸ obiektu oraz linie poziomów: pocz¹tkowego (y0) i ustalonego (y_ust)
p1 = plot(t, aTwew, 'b-', 'LineWidth', 2.2);
yline(y0, 'k--', 'Poziom pocz¹tkowy T_{wew0}', 'LabelHorizontalAlignment', 'left');
yline(y_ust, 'k--', 'Wartoœæ ustalona T_{wew,ust}', 'LabelHorizontalAlignment', 'left');
xline(czas_skok, 'r:', 'Chwila skoku t_{skok}', 'LabelVerticalAlignment', 'bottom');

% 2. Narysowanie linii stycznej (od przeciêcia z y0 do przeciêcia z y_ust)
t_styczna_vec = [t_start_st, t_koniec_st];
y_styczna_vec = [y0, y_ust];
p2 = plot(t_styczna_vec, y_styczna_vec, 'r-', 'LineWidth', 1.8);

% 3. Punkt przegiêcia oraz punkty przeciêcia stycznej z asymptotami
p3 = plot(t_przeg, y_przeg, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(t_start_st, y0, 'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(t_koniec_st, y_ust, 'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'k');

% 4. Rzutowanie punktu przeciêcia górnego na oœ czasu (do pokazania T)
plot([t_koniec_st, t_koniec_st], [y0, y_ust], 'k:', 'LineWidth', 1.2);

% 5. Podpisy przedzia³ów T0 oraz T na poziomie y0
text((czas_skok + t_start_st)/2, y0 + 0.08, sprintf('T_0 = %.1f s', T0_st), ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'r');
text((t_start_st + t_koniec_st)/2, y0 + 0.08, sprintf('T = %.1f s', T_st), ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', [0 0.5 0]);
text(t_przeg + 40, y_przeg, sprintf('Punkt przegiêcia\n(t = %.1f s, T = %.2f °C)', t_przeg, y_przeg), ...
     'FontSize', 9, 'FontWeight', 'bold', 'Color', 'r');

% 6. Krzywa wyznaczonego modelu ze stycznej
p4 = plot(t, Twew_model_st, 'g-.', 'LineWidth', 1.8);

ylim([y0 - 0.15, y_ust + 0.2]);
legend([p1, p2, p3, p4], {'OdpowiedŸ obiektu T_{wew}', 'Styczna w pkt. przegiêcia', ...
       'Punkt przegiêcia', 'Model z metody stycznej'}, 'Location', 'southeast');