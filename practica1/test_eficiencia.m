% Este script ejecuta los 3 métodos (Explícito, Implícito, C-N)
% para un rango de mallas (J_vals) y genera las gráficas 
% comparativas de eficiencia (Error vs. Tiempo).

% Genera dos gráficas: una para error relativo y otra para absoluto.

close all;

disp('--- INICIANDO TEST DE EFICIENCIA COMPARATIVA ---');

% --- 1. Configuración de la Simulación ---
J_vals = [100:100:500]; 
N_vals = J_vals * 10; 

error_types_names = {'Relativo', 'Absoluto'};

for tipo_error = 1:2
    
    nombre_error = error_types_names{tipo_error};

    fprintf('Calculando para: Error %s (choiceError = %d)\n', nombre_error, tipo_error);
    
    fprintf('  Ejecutando (1) Explícito...\n');
    % Para Explícito (1), pasamos [] a N_values para que use la cota interna
    [Times_Exp, Errors_Exp] = practica1_1(1, tipo_error, J_vals, []); 
    
    fprintf('  Ejecutando (2) Implícito...\n');
    [Times_Imp, Errors_Imp] = practica1_1(2, tipo_error, J_vals, N_vals);
    
    fprintf('  Ejecutando (3) Crank-Nicolson...\n');
    [Times_CN,  Errors_CN]  = practica1_1(3, tipo_error, J_vals, N_vals);

    figure_num = 3 + tipo_error; 
    figure(figure_num);
    clf;
    hold on;
    
    % La idea es graficar la eficiencia de cada metodo en colores diferentes para poder apreciar cual es el más eficiente
    % es decir, error vs tiempo en cada metodo, y asi porder ver cual es el mejor metodo, es decir aquel que este mas abajo 
    % a la izquierda de entre los métodos

    % Dibujar los puntos de cada método
    scatter(Times_Exp(:), Errors_Exp(:), 'r', 'filled', 'DisplayName', 'Explícito (Estable)');
    scatter(Times_Imp(:), Errors_Imp(:), 'b', 'filled', 'DisplayName', 'Implícito');
    scatter(Times_CN(:), Errors_CN(:), 'g', 'filled', 'DisplayName', 'Crank-Nicolson');
    
    hold off;
    
    set(gca, 'XScale', 'log', 'YScale', 'log'); % Escala log-log
    xlabel('Tiempo de cómputo (s)');
    ylabel(sprintf('Error máximo (%s)', nombre_error));
    title(sprintf('Comparación de Eficiencia P1 (Error %s)', nombre_error));
    legend('Location', 'best');
    grid off; % esto hace que se vea mejor
    
    % Crea un nombre de archivo dinámico
    filename = sprintf('P1_comparacion_eficiencia_GENERAL_%s.png', upper(nombre_error));
    print(sprintf('-f%d', figure_num), filename, '-dpng');
    
    fprintf('Gráfica generada: %s\n', filename);
end
