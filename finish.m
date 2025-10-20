% Automatisches Commit + Push beim Beenden von MATLAB

% Git-Verzeichnis
repoPath = 'C:\Users\HIN1EH\MATLAB-Backups';
% Wechsle ins Git-Verzeichnis
cd(repoPath); 

% Prüfen, ob Änderungen vorhanden sind
[~, cmdout] = system('git status --porcelain');

if ~isempty(cmdout)
    % Füge alle Dateien außer in .gitignore hinzu
    system('git add .');
    
    % Commit + Kommentar
    commitMessage = ['[AUTO] Auto-Commit beim Beenden von MATLAB ' datestr(now,'[yyyy-mm-dd; HH:MM:SS]')];
    system(['git commit -m "' commitMessage '"']);
    
    % Push zum Remote Repository
    system('git push origin master');
end
