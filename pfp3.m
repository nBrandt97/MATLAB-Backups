clear; clc; close all;

% Wahrer Sollverlauf
x_true = linspace(0, 20, 200);
y_true = 0.5*x_true.^3 - 2*x_true.^2 + 3*x_true + 5;

% Stützstellen
x = linspace(0, 20, 5);  
y_meas = 0.5*x.^3 - 2*x.^2 + 3*x + 5 + randn(size(x))*200;  % Messpunkte mit Rauschen

% Startpolynom (Ordnung = Ordnung Soll?)
p0 = [0.5, -2, 3, 5];  

% Ordnung des Residuen-Fits
n_fit = 3;  

% Startpolynom auswerten
y_p0 = polyval(p0, x_true);

% Residuen berechnen
y_fit = y_meas - polyval(p0, x);    

% Fit der Residuen
p_fit = polyfit(x, y_fit, n_fit);

% Residuen-Fit auswerten
y_res_fit = polyval(p_fit, x_true);

% polyfitplus = Startpolynom + Residuen-Fit
y_plus = y_p0 + y_res_fit;

% Klassischer polyfit für Vergleich
n_total = 3;
p_classic = polyfit(x, y_meas, n_total);
y_classic = polyval(p_classic, x_true);

% Plot
figure; hold on; grid on;
plot(x_true, y_true, 'm-', 'LineWidth', 2);     % Sollverlauf
plot(x, y_meas, 'ro', 'MarkerFaceColor', 'r');  % Messpunkte
plot(x_true, y_plus, 'k-', 'LineWidth', 2);     % polyfitplus
plot(x_true, y_classic, 'g--', 'LineWidth', 2); % klassischer polyfit
plot(x_true, y_p0, 'b-.', 'LineWidth', 1.5);    % Startpolynom

legend('Sollwert (wahr)', 'Messpunkte', 'polyfitplus', 'klassischer polyfit', 'Startpolynom p0', 'Location', 'NorthWest');
xlabel('x'); ylabel('y');
title('polyfitplus vs. klassischer polyfit');

