clear; clc; close all;

% -------------------
% Parameter
% -------------------
a_mean = 1; a_std = 0.05; b = 10; c = 5;

% Monte-Carlo-Durchläufe
N = input('Wie viele Monte-Carlo-Durchläufe möchten Sie durchführen? ');
if isempty(N), N = 100; end
fprintf('Anzahl der Monte-Carlo-Durchläufe: %d\n', N);

% Anzahl Stützstellen für Polyfit (mind. 3)
supportPointsList = [3,4,5];

% Figur
figure; hold on; grid on;
xlabel('x'); ylabel('y');
title('Monte-Carlo-Parabeln und Polyfits');

% Daten erzeugen
[x, y, a_true, hMC] = plotMonteCarloParabeln(N, a_mean, a_std, b, c);

% Polyfits plotten
hPF = plotPolyfits(x, y, N, supportPointsList);

legend([hMC(1), hPF(1)], {'Monte-Carlo-Parabeln', 'Polyfit-Beispiel'}, 'Location','best');

% -------------------
% Funktionen
% -------------------

function [x, y, a, hMC] = plotMonteCarloParabeln(N, a_mean, a_std, b, c)
    x = linspace(-10, 10, 20);
    a = a_mean + a_std*randn(N,1);
    y = zeros(N, length(x));
    for k = 1:N
        y(k,:) = a(k)*x.^2 + b*x + c;
    end
    hMC = plot(x, y', 'b'); % jede Zeile -> eine Kurve
end

function hPF = plotPolyfits(x, y, N, supportPointsList)
    colors = {'r--','g--','m--'};
    hPF = gobjects(numel(supportPointsList),1);
    for jIdx = 1:numel(supportPointsList)
        j = supportPointsList(jIdx);
        idx = round(linspace(1, length(x), j));
        x_sub = x(idx);
        % Nur einige wenige Fits zeichnen, um Übersicht zu behalten
        for k = 1:min(N,10)
            y_sub = y(k, idx);
            [p, S, mu] = polyfit(x_sub, y_sub, 2);
            x_fit = linspace(min(x_sub), max(x_sub), 100);
            y_fit = polyval(p, x_fit, S, mu);
            hPF(jIdx) = plot(x_fit, y_fit, colors{jIdx});
        end
    end
end
