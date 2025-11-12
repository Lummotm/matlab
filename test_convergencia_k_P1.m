% Test de convergencia K (TIEMPO) - PARTE 1 (u_t = u_xx)
close all;
clear;
disp('--- INICIANDO TEST CONVERGENCIA K (Parte 1) ---');

T = 0.5;
J_values = [2000]; 
N_values = [100, 200, 400, 800];
k_values = T ./ N_values; 
method_names = {'Explicito', 'Implicito', 'Crank-Nicolson'};

% No trabajamos con el explicito, ya que el ese metodo fuerza la necesidad del ajuste de mu < 0.5 
% lo que generaria datos falsos al estar trabajando E~O(k+h^2) y con k=0.4*h^2 E~O(0.4*h^2+h^2)~O(h^2)
% lo que no verifica el orden de convergencia en tiempo sino el de espacio que ya conociamos antes 
for i = 2:3 
    choiceMethod = i; 
    choiceError = 1; 
    
    fprintf("P1 (k): Ejecutando método %s...\n", method_names{i});

    [~, Errors, ~, ~, ~] = practica1_1(choiceMethod, choiceError, J_values, N_values);

    % Extracción de datos (J está fijo, k varía)
    Errors_k = Errors(:, 1);
    
    % Ajuste por polyfit
    p_k = polyfit(log(k_values), log(Errors_k), 1);
    orden_q = p_k(1);
    fprintf('El orden de convergencia temporal (q) es: %.4f\n', orden_q);

    
    figure;

    % Graficamos nuestros datos
    loglog(k_values, Errors_k, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;

    % Graficamos la línea de ajuste
    k_fit = logspace(log10(min(k_values)), log10(max(k_values)), 20);
    Error_fit = exp(p_k(2)) * k_fit.^p_k(1);
    loglog(k_fit, Error_fit, 'r--', 'LineWidth', 1.5);

    % Leyenda
    title(sprintf('P1 Convergencia Temporal (q) - %s', method_names{i}));
    xlabel('log(k)'); 
    ylabel('log(Error)');
    legend('Datos', sprintf('Ajuste (q = %.4f)', orden_q), 'Location', 'best');
    grid on;

    print(sprintf('P1_convergencia_k_Metodo_%d.png', i), '-dpng');
end
