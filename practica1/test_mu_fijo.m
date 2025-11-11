
% Test de convergencia para verificar convergencia de h

close all;

disp('--- INICIANDO TEST DE CONVERGENCIA DE H---');

for i = 1:3 
    choiceMethod = i; 
    choiceError = 1; 
    J_values = [200:200:1000]; 
    if i == 1
        N_values = [];
    else
        N_values = J_values; 
    end
    h_values = 1./J_values;
    fprintf("Ejecutando para método %i\n",i)

    [Times, Errors, ~, ~, ~] = practica1_1(choiceMethod, choiceError, J_values, N_values);

    % Extraemos SOLO la diagonal, que tiene mu=0.4 constante
    Errors_diag = diag(Errors);
    Times_diag = diag(Times); 

    % Hacemos polyfit entre log(h) y log(Error) de la diagonal
    p = polyfit(log(h_values), log(Errors_diag), 1);

    orden_convergencia = p(1);

    fprintf('El orden de convergencia en espacio (p) es: %.4f\n', orden_convergencia);


    % % Graficamos
    figure;
    loglog(h_values, Errors_diag, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;

    % Dibuja la línea de ajuste
    h_fit = logspace(log10(min(h_values)), log10(max(h_values)), 20);
    Error_fit = exp(p(2)) * h_fit.^p(1);
    loglog(h_fit, Error_fit, 'r--', 'LineWidth', 1.5);

    title(sprintf('Convergencia Espacial (p) - Método %d', i));
    xlabel('log(h)');
    ylabel('log(Error)');
    legend('Datos', sprintf('Ajuste (p = %.4f)', orden_convergencia), 'Location', 'best');
    grid on;

    % print(sprintf('P1_convergencia_h_Metodo_%d.png', i), '-dpng');

end
