clear; clc; close all;

% Parameter
a_mean = 1; a_std = 0.005; b = 10; c = 5;

% Benutzerabfrage
N = input('Wie viele Monte-Carlo-Durchläufe möchten Sie durchführen? ');
if isempty(N), N = 100; end
fprintf('Anzahl der Monte-Carlo-Durchläufe: %d\n', N);

% Polyfit-Stützpunkte
% Definiert eine Liste mit verschiedenen Anzahlen an Stützstellen (zur Untersuchung des Fits).
supportPointsList = [2,3,4,5];

% Erzeuge neue Figur
figure;

% Funktionsaufruf
[x, y, hMC] = plotMonteCarloParabeln(N, a_mean, a_std, b, c);
hPF = plotPolyfits(x, y, N, supportPointsList);

% Übergebe Handles an Legende für korrekte Linienart in Legende
legend([hMC(1), hPF], {'Monte-Carlo-Parabeln', 'Polyfits'});

% Funktionen
% Funktion für Monte-Carlo-Parabeln
function [x, y, hMC] = plotMonteCarloParabeln(N, a_mean, a_std, b, c)
    x = linspace(-10, 10, 20);
    a = a_mean + a_std*randn(N,1);
    y = a .* (x.^2) + b*x + c;
    
    hMC = plot(x, y, 'b');   % Plot und Handle für Legende
    hold on;
    grid on;
    xlabel('x'); ylabel('y'); title('Funktionsschar');
end

% Funktion für Polyfits
function hPF = plotPolyfits(x, y, N, supportPointsList)
    for j = supportPointsList
        idx = round(linspace(1, length(x), j));
        x_subset = x(idx);
       
        fprintf('\nPolyfit mit %d Stützstellen:\n', j);
        for k = 1:N
            y_subset = y(k, idx);
            p = polyfit(x_subset, y_subset, 2);
            fprintf('  Durchlauf %d: a_fit = %.4f, b_fit = %.0f, c_fit = %.0f\n', ...
                    k, p(1), p(2), p(3));
            
            x_fit = linspace(min(x_subset), max(x_subset), 50);
            y_fit = polyval(p, x_fit);
            
            hPF = plot(x_fit, y_fit, 'r--'); % Handle für Legende wird überschrieben, reicht
        end
    end
end

%Test
%Test2
%Test3