clear; clc; close all;

% -------------------
% Parameter
% -------------------
a_mean = 1; a_std = 0.05; b = 10; c = 5;

% Benutzerabfrage
N = input('Wie viele Monte-Carlo-Durchläufe möchten Sie durchführen? ');
if isempty(N) 
    N = 100; 
end
fprintf('Anzahl der Monte-Carlo-Durchläufe: %d\n', N);

% Polyfit-Stützpunkte (Anzahlen)
supportPointsList = [2,3,4,5];

% Erzeuge Figur
figure('Units','normalized','Position',[0.05 0.05 0.9 0.8]);

% Funktionsaufruf: Monte-Carlo-Parabeln (liefert x, y und a_true)
[x, y, a_true, hMC] = plotMonteCarloParabeln(N, a_mean, a_std, b, c);

% Polyfits (liefert Handles)
hPF = plotPolyfits(x, y, a_true, N, supportPointsList);

% Legende: ein Beispiel-Handle aus MonteCarlo und ein Handle aus Polyfits
% (hPF ist Zelle mit Handles pro j; nehme das erste nicht-empty)
pfHandle = [];
for i = 1:numel(hPF)
    if ~isempty(hPF{i})
        pfHandle = hPF{i}(1);
        break;
    end
end
legend([hMC(1), pfHandle], {'Monte-Carlo-Parabeln', 'Polyfits (Beispiel)'}, 'Location','best');

% -------------------
% Funktionen
% -------------------

function [x, y, a, hMC] = plotMonteCarloParabeln(N, a_mean, a_std, b, c)
    % x: 1 x M
    x = linspace(-10, 10, 20);      % fein genug; M = 20
    M = length(x);
    % a: N x 1 (echte Koeffizienten)
    a = a_mean + a_std*randn(N,1);
    % y: N x M, jede Zeile eine Parabel
    y = zeros(N, M);
    for k = 1:N
        y(k, :) = a(k) * (x.^2) + b * x + c;
    end

    % Plot: MATLAB kann eine Matrix plotten: jede Zeile -> separate Kurve
    hMC = plot(x, y');   % transponiert: x ist 1xM, y' ist MxN -> jede Linie
    hold on;
    grid on;
    xlabel('x'); ylabel('y'); title('Monte-Carlo Parabeln und Polyfits');
end

function hPF = plotPolyfits(x, y, a_true, N, supportPointsList)
    % x: 1 x M
    % y: N x M
    M = length(x);
    colors = {'r--','g--','m--','c--'}; % nur Muster (wird zyklisch verwendet)
    hPF = cell(size(supportPointsList));
    jIndex = 0;
    for j = supportPointsList
        jIndex = jIndex + 1;
        handles_this_j = gobjects(0);
        % Erzeuge Indizes (Achtung: runden kann in extremen Fällen doppelte Indices erzeugen)
        idx = round(linspace(1, M, j));
        idx = max(1, min(M, idx)); % Sicherheit
        if numel(unique(idx)) < numel(idx)
            warning('Beim j=%d wurden doppelte Indices erzeugt — überprüfe idx: %s', j, mat2str(idx));
            idx = unique(idx,'stable');
        end
        x_subset = x(idx);  % 1 x j (geordnet)
        fprintf('\nPolyfit mit %d Stützstellen (Indices: %s):\n', j, mat2str(idx));

        % wir plotten nicht alle N Fits, sondern maximal maxPlotsPerJ (sonst Überlagerung)
        maxPlotsPerJ = min(20, N); % du kannst anpassen
        % wähle zufällige Teilmenge, damit auch bei großem N etwas Varianz sichtbar ist
        plotIdxSample = randperm(N, min(N, maxPlotsPerJ));

        for k = 1:N
            y_subset = y(k, idx); % 1 x j

            % Numerisch stabilere polyfit-Variante: gib mu zurück (centering+scaling)
            % p(1) entspricht quadratischem Koeffizienten (a_fit)
            [p, S, mu] = polyfit(x_subset, y_subset, 2);

            % Debug-Ausgabe des Fits (nur grob)
            a_fit = p(1);
            b_fit = p(2);
            c_fit = p(3);

            % Falls a_true übergeben, zeige Abweichung; ansonsten skip
            if nargin >= 3 && ~isempty(a_true)
                a_err = a_fit - a_true(k);
            else
                a_err = NaN;
            end
            % Zeige nur für die ersten paar Durchläufe detailliert, ansonsten kürzer
            if k <= 5
                fprintf('  Durchlauf %3d: a_fit = %.5f (err=%.5f), b_fit = %.3f, c_fit=%.3f\n', ...
                        k, a_fit, a_err, b_fit, c_fit);
            elseif k==6
                fprintf('  ... (weitere Durchläufe werden nicht vollständig geloggt)\n');
            end

            % Erzeuge feine x-Achse für Fit-Kurve
            x_fit = linspace(min(x_subset), max(x_subset), 80);
            y_fit = polyval(p, x_fit, S, mu); % numerisch stabil ausgewertet

            % Plot nur eine Stichprobe von Fits zur Übersicht
            if ismember(k, plotIdxSample)
                h = plot(x_fit, y_fit, colors{mod(jIndex-1,numel(colors))+1});
                handles_this_j(end+1) = h; %#ok<AGROW>

                % Plotte die Stützpunkte als Marker für diesen Durchlauf (klein)
                plot(x_subset, y_subset, 'ko', 'MarkerFaceColor','k', 'MarkerSize',4);
            end

            % ----- Optional: automatischer Qualitätscheck -----
            % Falls a_err auffällig groß ist, markiere und gib Details aus
            if ~isnan(a_err) && abs(a_err) > 0.05 % Schwellenwert anpassbar
                % zeichne die abhängigen Werte fett (nur für die sample-Kurven)
                if ismember(k, plotIdxSample)
                    plot(x_fit, y_fit - 0.02*abs(y_fit), '--', 'LineWidth', 1.2); % visuelle Markierung
                end
                fprintf('    --> WARNUNG: Durchlauf %d hat großen a-Fehler: %.5f\n', k, a_err);
            end
        end

        hPF{jIndex} = handles_this_j;
    end
end
