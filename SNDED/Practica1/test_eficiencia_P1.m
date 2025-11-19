% Este test genera la gráfica de eficiencia con todos los nodos de (error,tiempo) de cada método.

close all;
clear;
disp(' INICIANDO TEST DE EFICIENCIA (P1 COMPLETO) ');

% Configuración 
% Valores de N y J para Implícito y C-N
J_vals = [50, 100, 150, 200, 250]; 
N_vals = [100, 500, 1000, 2000, 4000, 8000]; 

% J para el Explícito
J_vals_exp = [50, 100, 150, 200, 250]; 

[Times_Exp, Errors_Exp] = practica1_1(1, J_vals_exp, []); 

[Times_Imp, Errors_Imp] = practica1_1(2, J_vals, N_vals);

[Times_CN,  Errors_CN]  = practica1_1(3, J_vals, N_vals);

% Creación de la Gráfica Principal 
figure_num = 10; % Un numero de figura que no se solape con las otras
figure(figure_num);
clf;
hold on;

scatter(Times_Exp(:), Errors_Exp(:), 'r', 'filled', 'DisplayName', 'Explícito (Estable)');
scatter(Times_Imp(:), Errors_Imp(:), 'b', 'filled', 'DisplayName', 'Implícito');
scatter(Times_CN(:), Errors_CN(:), 'g', 'filled', 'DisplayName', 'Crank-Nicolson');

hold off;
set(gca, 'XScale', 'log', 'YScale', 'log'); % Escala log-log
xlabel('Tiempo de cómputo (s)');
ylabel('Error máximo (Relativo)');
title('Comparación de Eficiencia P1');
legend('Location', 'best');

print(sprintf('-f%d', figure_num), 'P1_comparacion_eficiencia_GENERAL_RELATIVO.png', '-dpng');

