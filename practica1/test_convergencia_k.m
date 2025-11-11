% Test de convergencia para verificar convergencia de k

close all;
T = 0.5;

disp('--- INICIANDO TEST DE CONVERGENCIA DE K---');

for i = 2:3 
    choiceMethod = i; 
    choiceError = 1; 
    J_values = [2000]; % J muy grande para aislar el error temporal
    N_values = [100, 200, 400, 800];
    k_values = T ./ N_values; 
    fprintf("Ejecutando para método %i\n",i)

    [Times, Errors, ~, ~, ~] = practica1_1(choiceMethod, choiceError, J_values, N_values);


    % Hacemos polyfit entre log(k) y log(Error)
    % p_k(1) será la pendiente (el orden q)
    p_k = polyfit(log(k_values), log(Errors(:, 1)), 1);
    orden_q = p_k(1);
    fprintf('El orden de convergencia temporal (q) es: %.4f\n', orden_q);

    % figure;
    % % Ploteamos k_values vs Errors(:, 1)
    % loglog(k_values, Errors(:, 1), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    % hold on;
    %
    % % Dibuja la línea de ajuste
    % k_fit = logspace(log10(min(k_values)), log10(max(k_values)), 20);
    % Error_fit = exp(p_k(2)) * k_fit.^p_k(1);
    % loglog(k_fit, Error_fit, 'r--', 'LineWidth', 1.5);
    %
    % title(sprintf('Convergencia Temporal (q) - Método %d', i));
    % xlabel('log(k)'); 
    % ylabel('log(Error)');
    % % Variable de leyenda corregida (orden_q)
    % legend('Datos', sprintf('Ajuste (q = %.4f)', orden_q), 'Location', 'best');
    % grid on;
    %
    % % Guardar la gráfica
    % print(sprintf('P1_convergencia_k_Metodo_%d.png', i), '-dpng');
end
