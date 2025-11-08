% Este script ejecuta los 3 métodos (Explícito, Implícito, C-N)
% para un rango de mallas (J_vals) y genera las gráficas 
% comparativas de eficiencia (Error vs. Tiempo).

% Genera dos gráficas: una para error relativo y otra para absoluto.

clc;
clear all;
close all;

disp('--- INICIANDO TEST DE EFICIENCIA COMPARATIVA ---');

% --- 1. Configuración de la Simulación ---
J_vals = [100:100:500]; 
N_gen = J_vals * 10; 

error_types_names = {'Relativo', 'Absoluto'};

for tipo_error = 1:2
    
    nombre_error = error_types_names{tipo_error};

    fprintf('Calculando para: Error %s (choiceError = %d)\n', nombre_error, tipo_error);
    
    fprintf('  Ejecutando (1) Explícito...\n');
    % Para Explícito (1), pasamos [] a N_values para que use la cota interna
    [Times_Exp, Errors_Exp] = practica1_1(1, tipo_error, J_vals, []); 
    
    fprintf('  Ejecutando (2) Implícito...\n');
    [Times_Imp, Errors_Imp] = practica1_1(2, tipo_error, J_vals, N_gen);
    
    fprintf('  Ejecutando (3) Crank-Nicolson...\n');
    [Times_CN,  Errors_CN]  = practica1_1(3, tipo_error, J_vals, N_gen);

    figure_num = 3 + tipo_error; 
    figure(figure_num);
    clf;
    hold on;
    
    % Dibujar los puntos de cada método
    scatter(Times_Exp(:), Errors_Exp(:), 'r', 'filled', 'DisplayName', 'Explícito (Estable)');
    scatter(Times_Imp(:), Errors_Imp(:), 'b', 'filled', 'DisplayName', 'Implícito');
    scatter(Times_CN(:), Errors_CN(:), 'g', 'filled', 'DisplayName', 'Crank-Nicolson');
    
    hold off;
    
    set(gca, 'XScale', 'log', 'YScale', 'log'); % Escala log-log
    xlabel('Tiempo de cómputo (s)');
    ylabel(sprintf('Error máximo (%s)', nombre_error));
    title(sprintf('Comparación de Eficiencia (Error %s)', nombre_error));
    legend('Location', 'best');
    grid off; % esto hace que se vea mejor
    
    % Crea un nombre de archivo dinámico
    filename = sprintf('comparacion_eficiencia_GENERAL_%s.png', upper(nombre_error));
    print(sprintf('-f%d', figure_num), filename, '-dpng');
    
    fprintf('-> Gráfica generada: %s\n', filename);

end
