clear; clc; close all;

% gegebene Funktion: y = ax^2 + bx + c

% Deklaration & Initialisierung Variablen
a_mean = 1;
a_std = 0.005;
b = 10;
c = 5;

% Benutzerabfrage
N = input('Wie viele Monte-Carlo-Durchläufe möchten Sie durchführen? ');
if isempty(N)                          % falls der Benutzer nur Enter drückt
    N = 100;                           % Standardwert
end

fprintf('Anzahl der Monte-Carlo-Durchläufe: %d\n', N);

% erzeuge 20 gleichmäßige Stützpunkte zwischen -10 & 10
x = linspace(-10, 10, 20);             

% erzeuge Spaltenvektor mit normalverteilten Zufallswerten + Anpassung an
% gegegebene Bedingungen
a = a_mean + a_std * randn(N,1);       

figure;                                % neues Plotfenster
hold on;    
grid on;
xlabel('x'); 
ylabel('y');
title('Funktionsschar');% behalte alte Plots bei neuer Schleifeniteration

% Schleife für Funktionsschar
for i = 1:N                            % Schleife von 1 bis N
y = a(i)*x.^2 + b*x + c;               % berechne y mit i-ten Eintrag von a
plot(x,y);                             % erzeuge Plot
end

% Ploteigenschafen
% Anzahl der Stützstellen
supportPointsList = [2,3,4,5];

for j = supportPointsList
    fprintf('Polyfit mit %d Stützstellen:\n', j);

    % zufällige Auswahl von j Stützstellen aus x
    idx = round(linspace(1, length(x), j));     % erzeuge Vektor mit ganzzahligen gleichmäßigen j-Einträgen zw. 1 bis 20
    x_subset = x(idx);                          % wähle aus Vektor x Komponenten die idx entsprechen

    % Schleife über Monte-Carlo-Durchläufe
    for k = 1:N
        y = a(k)*x.^2 + b*x + c;

        y_subset = y(idx);                      % y-Werte für die ausgewählten Stützstellen
        p = polyfit(x_subset, y_subset, 2);     % Quadratischer Fit

        % Ausgabe gefittete Koeffizienten
        fprintf('  Durchlauf %d: a_fit = %.4f, b_fit = %0.0f, c_fit = %0.0f\n', ... 
            k, p(1), p(2), p(3));

        
        x_fit = linspace(min(x_subset), max(x_subset), 100);
        y_fit = polyval(p, x_fit);
        plot(x_fit, y_fit, '--');               % gestrichelter Fit
    end
end


