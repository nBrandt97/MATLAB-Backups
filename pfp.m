%% polyfit vs. polyfitplus Beispiel mit Sollwertverlauf

clear; clc; close all;

%% 1. Sollwertverlauf (ideal)
x_dense = linspace(0, 5, 200);          % dichte x-Werte für Kurve
y_soll = 0.05*x_dense.^5 - 0.2*x_dense.^4 + x_dense.^3 - 0.5*x_dense.^2 + 1.5*x_dense + 2;

%% 2. Messpunkte (z.B. simuliert mit kleinen Abweichungen)
x_meas = [0.5 1.5 2.5 3 3.5];          % Stützstellen
y_meas = interp1(x_dense, y_soll, x_meas) + 0.2*randn(size(x_meas));  % kleine Messabweichungen

%% 3. Startpolynom p0 (z.B. theoretisch aus Teilinformationen, Grad 4)
p0 = [0.1, -0.2, 1, 0.5, 1]; % Grad 4: a4..a0

%% 4. polyfitplus
n = 5;       % gewünschter Endgrad
nfit = 2;    % nur Terme bis 2. Ordnung fitbar

% Startpolynom auswerten an Messpunkten
yp0 = polyval(p0, x_meas);

% Differenz zur Messung bilden
yfit = y_meas - yp0;

% Teilfit für niedrige Terme
pfit = polyfit(x_meas, yfit, nfit);

% pfit auffüllen auf Länge p0
pfit_extended = [zeros(1, length(p0)-length(pfit)), pfit];

% Endpolynom
p_plus = p0 + pfit_extended;

%% 5. Klassischer polyfit (Gesamtfit)
nklass = min(n, length(x_meas)-1);  % klassische polyfit Einschränkung
p_classic = polyfit(x_meas, y_meas, nklass);

%% 6. Auswertung für dichten Plot
y_plus = polyval(p_plus, x_dense);
y_classic = polyval(p_classic, x_dense);

%% 7. Plotten
figure; hold on; grid on;
plot(x_dense, y_soll, 'm-', 'LineWidth', 2);       % Sollwertverlauf
plot(x_meas, y_meas, 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r'); % Messpunkte
plot(x_dense, y_plus, 'k-', 'LineWidth', 2);       % polyfitplus
plot(x_dense, y_classic, 'g-.', 'LineWidth', 2);  % klassischer polyfit
plot(x_dense, polyval(p0, x_dense), 'b--', 'LineWidth', 1.5);         % Startpolynom

legend('Sollwertverlauf','Messpunkte','polyfitplus','klassischer polyfit','Startpolynom p0','Location','NorthWest');
xlabel('x'); ylabel('y');
title('polyfitplus vs. klassischer polyfit mit Sollwertverlauf');
