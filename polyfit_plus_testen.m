    clear
     %Daten für Sensor aus Spalte 1
     %ID=1.903174937577550e+14
     load('\\bosch.com\dfsrb\DfsDE\LOC\Eh\GS\GS-SI_ENS-Eh\GS_EZS-Eh\Ma\Praktikanten\Praktikant Feldvoß\02_HPS5_Abgleich_hörer_Ordnung\03_Matlab\polyfit_plus_test\druck_lin_test_.mat')
     load('\\bosch.com\dfsrb\DfsDE\LOC\Eh\GS\GS-SI_ENS-Eh\GS_EZS-Eh\Ma\Praktikanten\Praktikant Feldvoß\02_HPS5_Abgleich_hörer_Ordnung\03_Matlab\polyfit_plus_test\p_out_calib.mat')
                  degree1 = 3  
     nonlinfit_plus = polyfit_plus(druck_lin_test_,p_out_calib,degree1,[2.466729946867246e-04,2.757336357207350e-05];
                  %Reihenfolge von p0=[ CLin2_mittel CLin3_mittel] Achtung automatisch in dieser Funktion hier übergeben!
                  %nonlinfit_plus= [höchster Grad.....niedrigster Grad]
                  CLin0_C_undis = nonlinfit_plus(degree1+1);
                  CLin1_C_undis = nonlinfit_plus(degree1);
                  CLin2_C_undis = nonlinfit_plus(degree1-1);
                  CLin3_C_undis = nonlinfit_plus(degree1-2); 