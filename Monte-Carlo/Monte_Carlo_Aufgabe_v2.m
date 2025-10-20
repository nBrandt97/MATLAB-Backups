clear; clc; close all;

% gegebene Funktion: y = ax^2 + bx + c

% Parameter
a_mean = 1; a_std = 0.005; b = 10; c = 5;

% Benutzerabfrage
N = input('Wie viele Monte-Carlo-Durchläufe möchten Sie durchführen? ');
if isempty(N) 
    N = 100; 
end
fprintf('Anzahl der Monte-Carlo-Durchläufe: %d\n', N);

% Stützpunkte
% Erzeugt 20 gleichmäßig verteilte x-Werte zwischen -10 und 10.
x = linspace(-10, 10, 20);  

% Monte-Carlo-Koeffizient
% Erstellt N zufällige Werte für a, normalverteilt um a_mean mit Standardabweichung a_std.
a = a_mean + a_std * randn(N,1);  

% Berechnung aller Parabeln
% Berechnet für jeden Zufallswert von a die y-Werte der Parabel an den 20 x-Stützpunkten.
% Ergebnis: eine Matrix (N x Länge(x)) mit allen Parabeln.
y = a .* (x.^2) + b * x + c;

% Plot Monte-Carlo-Parabeln
% Öffnet eine neue Grafik, behält alle Plots offen ("hold on"),
figure; 
hold on; 
grid on; 
xlabel('x');
ylabel('y');
plot(x, y);
title('Funktionsschar');

% Speichere ein Handle der Matrix der Parabeln (für Legende)
hMC = plot(x, y, 'b');

% Polyfit-Stützpunkte
% Definiert eine Liste mit verschiedenen Anzahlen an Stützstellen (zur Untersuchung des Fits).
supportPointsList = [2,3,4,5]; 

% Polyfit
for j = supportPointsList                                    % Schleife über verschiedene Anzahlen von Stützstellen.
    idx = round(linspace(1, length(x), j));                  % Wählt j gleichmäßig verteilte Indizes aus den x-Stützpunkten.
    x_subset = x(idx);                                       % Entnimmt die x-Werte an diesen Stützstellen.
   
    fprintf('\n Polyfit mit %d Stützstellen:\n', j);
    for k = 1:N                                              % Schleife über alle Monte-Carlo-Durchläufe.
        y_subset = y(k, idx);                                % Entnimmt die passenden y-Werte der aktuellen Parabel an den j Stützstellen.
        p = polyfit(x_subset, y_subset, 2);                  % Führt eine quadratische Polynomapproximation (Grad 2) durch.
        
        % Ausgabe der geschätzten Koeffizienten
        fprintf('  Durchlauf %d: a_fit = %.4f, b_fit = %.0f, c_fit = %.0f\n', ...
                k, p(1), p(2), p(3));
        
        % Plot Fit
        x_fit = linspace(min(x_subset), max(x_subset), 50);  % Erzeugt fein aufgelöste x-Werte im Bereich der Stützstellen.
        y_fit = polyval(p, x_fit);                           % Berechnet die y-Werte des Fits für diese x-Werte.
        hPF = plot(x_fit, y_fit, 'r--');                     % Speichere ein Handle für den letzten Polyfit
    end 
end

legend([hMC(1) hPF], {'Monte-Carlo-Parabeln', 'Polyfits'}); % Übergebe Handles an Legende für korrekte Linienart in Legende