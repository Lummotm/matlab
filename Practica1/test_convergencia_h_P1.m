% Test de convergencia h (ESPACIO) - PARTE 1 (u_t = u_xx)
close all;
clear;
disp(' INICIANDO TEST CONVERGENCIA H (Parte 1) ');

J_values = [100, 200, 300, 400, 500]; 
h_values = 1./J_values;
method_names = {'Explicito', 'Implicito', 'Crank-Nicolson'};

for i = 1:3 
    choiceMethod = i; 
    
    % Usamos el mu del ajuste automatico al pasar N vacio a mi función hace el ajuste de convergencia necesario
    N_values = []; 
    
    fprintf("P1 (h): Ejecutando método %s...\n", method_names{i});
    
    [~, Errors, ~, ~, ~] = practica1_1(choiceMethod, J_values, N_values);
    
    % Extraemos la diagonal, es decir forzamos mu = 0.1 cte.
    Errors_h = diag(Errors);
    
    % Ajuste por polyfit
    p = polyfit(log(h_values), log(Errors_h), 1);
    orden_p = p(1);
    fprintf('El orden de convergencia en espacio (p) es: %.4f\n', orden_p);
    
    figure; 
    
    % Graficamos datos
    loglog(h_values, Errors_h, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    
    % Graficamos ajuste via polyfit
    h_fit = logspace(log10(min(h_values)), log10(max(h_values)), 20);
    Error_fit = exp(p(2)) * h_fit.^p(1);
    loglog(h_fit, Error_fit, 'r--', 'LineWidth', 1.5);
    
    % Leyenda
    title(sprintf('P1 Convergencia Espacial (p) - %s', method_names{i}));
    xlabel('log(h)');
    ylabel('log(Error)');
    legend('Datos', sprintf('Ajuste (p = %.4f)', orden_p), 'Location', 'best');
    grid on;
    
    % Denotamos en función del método via sprintf
    print(sprintf('P1_convergencia_h_Metodo_%d.png', i), '-dpng');
end

%  Prueba Adicional: Implícito Sin Ajuste (P1) 
% Se repite para i=2 (Implícito) pero sin ajustar mu (N_values = J_values).
% Esto demuestra cómo el error temporal O(k) domina y reduce el orden p.
i = 2;
choiceMethod = i; 

% Usamos el mu del ajuste automatico al pasar N vacio a mi función hace el ajuste de convergencia necesario
N_values = J_values; 

fprintf("P1 (h): Ejecutando método %s...\n", method_names{i});

[~, Errors, ~, ~, ~] = practica1_1(choiceMethod, J_values, N_values);

% Extraemos la diagonal, es decir forzamos mu = 0.1 cte.
Errors_h = diag(Errors);

% Ajuste por polyfit
p = polyfit(log(h_values), log(Errors_h), 1);
orden_p = p(1);
fprintf('El orden de convergencia en espacio (p) sin ajustar es: %.4f\n', orden_p);

figure; 

% Graficamos datos
loglog(h_values, Errors_h, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;

% Graficamos ajuste via polyfit
h_fit = logspace(log10(min(h_values)), log10(max(h_values)), 20);
Error_fit = exp(p(2)) * h_fit.^p(1);
loglog(h_fit, Error_fit, 'r--', 'LineWidth', 1.5);

% Leyenda
title(sprintf('P1 Convergencia Espacial Sin Ajuste (p) - %s', method_names{i}));
xlabel('log(h)');
ylabel('log(Error)');
legend('Datos', sprintf('Ajuste (p = %.4f)', orden_p), 'Location', 'best');
grid on;

% Denotamos en función del método via sprintf
print(sprintf('P1_convergencia_h_Metodo_%d_no_ajustado.png', i), '-dpng');
