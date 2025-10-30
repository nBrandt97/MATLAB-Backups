clear; clc; close all;

% Wahrer Sollverlauf (5. Ordnung)
x_true = linspace(0, 20, 200);
y_true = 0.05*x_true.^5 - 0.2*x_true.^4 + x_true.^3 - 0.5*x_true.^2 + 1.5*x_true + 2;

% "Gemessene" Stützstellen – viel zu wenige!
x = linspace(0, 20, 6);  % nur 4 Punkte!
y_meas = 0.05*x.^5 - 0.2*x.^4 + x.^3 - 0.5*x.^2 + 1.5*x + 2 + randn(size(x))*50;

% Ziel: eigentlich 5. Ordnung, aber das geht mit polyfit nicht direkt
n_total = 5;     % gewünschte Gesamtordnung
nfit = 2;        % nur die unteren Terme werden gefittet

% Startpolynom (wir „wissen“ etwas über die höheren Terme)
p0 = [0.05, -0.2, 1, 0, 0, 0];  % Startpolynom bis Grad 5

% --- polyfitplus ---
yp0 = polyval(p0, x);
yfit = y_meas - yp0;
pfit = polyfit(x, yfit, nfit);
pfit_extended = [zeros(1, length(p0)-length(pfit)), pfit];
pplus = p0 + pfit_extended;

% --- klassischer polyfit (nur bis Grad 2 möglich) ---
p_classic = polyfit(x, y_meas, n_total);

% --- Auswertung ---
y_plus = polyval(pplus, x_true);
y_classic = polyval(p_classic, x_true);

% --- Plot ---
figure; hold on; grid on;
plot(x_true, y_true, 'm-', 'LineWidth', 2);            % Sollverlauf (wahr)
plot(x, y_meas, 'ro', 'MarkerFaceColor', 'r');         % Messpunkte
plot(x_true, y_plus, 'k-', 'LineWidth', 2);            % polyfitplus
plot(x_true, y_classic, 'g--', 'LineWidth', 2);        % klassischer polyfit
plot(x_true, polyval(p0, x_true), 'b-.', 'LineWidth', 1.5); % Startpolynom

legend('Sollwert (wahr)', 'Messpunkte', 'polyfitplus', ...
       'klassischer polyfit', 'Startpolynom p0', 'Location', 'NorthWest');
xlabel('x'); ylabel('y');
title('polyfitplus vs. klassischer polyfit – zu wenige Stützstellen');

