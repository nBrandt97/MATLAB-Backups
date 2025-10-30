clear; clc; close all;

% Parameter
n = 3;        % Gesamtgrad
nfit = 3;     % Grad des Teilfits

% Sollwertverlauf
x = linspace(0, 20, 100);
y = 0.05*x.^5 - 0.2*x.^4 + x.^3 - 0.5*x.^2 + 1.5*x + 2;

% Klassischer Polyfit
p = polyfit(x, y, n);

% Startpolynom (angenommen 4. Ordnung)
p0 = [0.1, -0.2, 1, 0.5, 1];

% Startpolynom auswerten
yp0 = polyval(p0, x);

% Differenz zu Sollwerten
yfit = y - yp0;

% Teilfit durchführen (nur bis nfit)
pfit = polyfit(x, yfit, nfit);

% Auf gleiche Länge bringen
pfit_extended = [zeros(1, length(p0)-length(pfit)), pfit];

% Endpolynom
pplus = p0 + pfit_extended;

% Werte für Plot
y_plus = polyval(pplus, x);

% Plot
figure; hold on; grid on;
plot(x, y, 'm-', 'LineWidth', 2);       % Sollwertverlauf
plot(x, yp0, 'b--', 'LineWidth', 1.5);  % Startpolynom
plot(x, y_plus, 'k-', 'LineWidth', 2);  % polyfitplus
plot(x, y, 'g-.', 'LineWidth', 2);      % klassischer polyfit

legend('Sollwertverlauf','Startpolynom p0','polyfitplus','klassischer polyfit','Location','NorthWest');
xlabel('x'); ylabel('y');
title('polyfitplus vs. klassischer polyfit mit Sollwertverlauf');
