function [Times, Errors, U_int, t, x] = practica1_1(choiceMethod, J_values, N_values)
    % El texto siguiente esta diseñado para que sea leido al usar help el esta función  
    % Function practica1_1. 
    % Resuelve la ecuación del calor 1D (u_t = u_xx)
    % Input: (choiceMethod, J_values, N_values)
    % - choiceMethod: metodo escogido 
    %       - 1) Explícito
    %       - 2) Implícito
    %       - 3) Crank-Nicolson
    % - J_values: valores de J, peudes pasar un vector y el programa se encarga de ello usar cada uno para cada malla de h
    % - N_values: valores de N, puedes pasar un escalar o un vector, si no se pasa nada, se considera el caso de mu = 0.4 < 1/2 para poder ejecutar el metodo explicito

    T = 0.5; % Tiempo final de la simulación

    % Condiciones iniciales y solución exacta 
    u0_fun = @(x) sin(2*pi*x);
    u_exact = @(x,t) sin(2*pi*x) .* exp(-4*pi^2 * t);

    if nargin < 1 || isempty(choiceMethod)
        choiceMethod = 0; % Valor por defecto si no se introduce
    end
    
    if nargin < 2 || isempty(J_values)
        J_values = [100]; 
    end

    if nargin < 3 || isempty(N_values)
        % Si no se da, se calcula N para CUMPLIR LA CONDICIÓN DE ESTABILIDAD
        % del método explícito (mu <= 1/2).
        % mu = k/h^2 = (1 / (2.5 * J^2)) / (1/J^2) = 1/2.5 = 0.4
        disp("INFO: 'N_values' no especificado. Calculando N por estabilidad (mu=0.4).")
        N_values = 2.5 * J_values.^2 * T;
    end

    % Bucle 'while' para asegurar que el usuario elige un método válido.
    validMethods = 1:3;
    while ~ismember(choiceMethod, validMethods)
        disp("")
        disp("Escoja un método.");
        disp("Métodos de resolución disponibles:");
        disp("  1) Explícito");
        disp("  2) Implícito");
        disp("  3) Crank-Nicolson");
        choiceMethod = input("Selecciona un método [1-3]: ");
        disp("")
    end

    % Se usa para los nombres de las graficas que se guardaran posteriormente
    method_names = {'Explícito', 'Implícito', 'Crank-Nicolson'};

    % Inicialización
    L_j = length(J_values);
    L_n = length(N_values); 
    h_values = 1 ./ J_values;
    k_values = T ./ N_values; 

    Times = zeros(L_n, L_j);
    Errors = zeros(L_n, L_j);
    U_int = [];

    % Iteramos
    for n = 1:L_n
        k = k_values(n);
        N = N_values(n);
        for j = 1:L_j
            J = J_values(j);
            h = h_values(j);
            mu = k / h^2;

            % Malla
            t = linspace(0, T, N+1)';
            x = linspace(0, 1, J+1)';
            x_int = x(2:end-1); 

            % Condición inicial
            u0 = u0_fun(x_int);

            % Meter lo parametros en un struct para facilitar lectura.
            params.J = J; 
            params.N = N; 
            params.mu = mu; 
            params.x = x_int; 
            params.t = t; 

            % Resolver
            [time, errors, U_int_current] = solver_selection(choiceMethod, params, u0, u_exact);

            errmax = max(errors);
            Times(n, j) = time;
            Errors(n, j) = errmax;

            % Guardamos la matriz, si se da un solo argumento para N,J
            if L_n == 1 && L_j == 1
                U_int = U_int_current;
            end
        end
    end

    nombre_metodo = method_names{choiceMethod};

    % Solo tiene caso iterar sobre casos de al menos 2 elementos para graficar, porque sino tendriamos elementos unipuntuales en las graficas que no aportan nada
    if  length(N_values) > 1 || length(J_values) > 1
        % 1) k fijo, variando h
        figure(1);
        clf;      
        hold on; 

        legend_text_k = cell(L_n, 1); 

        for n = 1:L_n
            loglog(Times(n,:), Errors(n,:), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
            legend_text_k{n} = sprintf('k = %.0e', k_values(n));
        end

        hold off; 
        grid on;
        xlabel('Tiempo de ejecución');
        ylabel('Error máximo');
        title(sprintf('(k fijo, var h): %s', nombre_metodo));
        legend(legend_text_k, 'Location', 'best');
        print("-f1", "P1_k_fijo_var_h_" + nombre_metodo, "-dpng");

        % 2) h fijo, variando k
        figure(2);
        clf;      
        hold on; 

        legend_text_h = cell(L_j, 1); 

        for j = 1:L_j
            loglog(Times(:,j), Errors(:,j), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
            legend_text_h{j} = sprintf('h = %.2e', h_values(j));
        end

        hold off; 
        grid on;
        xlabel('Tiempo de ejecución');
        ylabel('Error máximo');
        title(sprintf('(h fijo, var k): %s', nombre_metodo));
        legend(legend_text_h, 'Location', 'best');
        print("-f2", "P1_h_fijo_var_k_" + nombre_metodo, "-dpng");
    end
end


function [time, errors, U_int] = solver_selection(choiceMethod, params, u0, u_exact)
    switch choiceMethod
    case 1 % Explícito
        % Requiere condición: mu <= 1/2 para estabilidad
        if params.mu <= 1/2 
            tic;
            [U_int, errors] = solve_explicito(params, u0, u_exact);
            time = toc;
        else 
            % Inestable: devolvemos NaN
            errors(:) = NaN;
            time = NaN;
            U_int = NaN;
        end
    case 2 % Implícito 
        tic;
        [U_int, errors] = solve_implicito(params, u0, u_exact);
        time = toc;
    case 3 % Crank-Nicolson 
        tic;
        [U_int, errors] = solve_crank(params, u0, u_exact);
        time = toc;
    end
end

function [U_int, errors] = solve_explicito(params, u, u_exact)
    J = params.J;
    N = params.N;
    mu = params.mu;
    x = params.x;
    t = params.t;

    A = spdiags([mu , (1-2*mu), mu], -1:1, J-1, J-1);

    % Inicialización de la matriz de solución U (puntos interiores)
    U_int = zeros(J-1, N+1);
    U_int(:, 1) = u; 

    % Inicialización del error
    errors = zeros(1, N+1);
    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex)) / (max(u_ex) + eps) ;

    for n = 1:N
        u = A * u;
        U_int(:, n+1) = u; % Almacenar el paso de tiempo
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex)) / (max(u_ex) + eps);
    end
end

function [U_int, errors] = solve_implicito(params, u, u_exact)
    J = params.J;
    N = params.N;
    mu = params.mu;
    x = params.x;
    t = params.t;

    A = spdiags([-mu , (1+2*mu), -mu ], -1:1, J-1, J-1);
    dA = decomposition(A);

    % Inicialización de la matriz de solución U (puntos interiores)
    U_int = zeros(J-1, N+1);
    U_int(:, 1) = u; 

    % Inicialización del error
    errors = zeros(1, N+1);
    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex)) / (max(u_ex) + eps) ;

    for n = 1:N
        u = dA \ u;
        U_int(:, n+1) = u; 
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex)) / (max(u_ex) + eps);
    end
end

function [U_int, errors] = solve_crank(params, u, u_exact)
    J = params.J;
    N = params.N;
    mu = params.mu;
    x = params.x;
    t = params.t;

    temp = mu / 2;
    A = spdiags([-temp*ones(J-1,1), (1+mu)*ones(J-1,1) , -temp*ones(J-1,1) ], -1:1, J-1, J-1);
    B = spdiags([temp*ones(J-1,1),  (1-mu)*ones(J-1,1),  temp*ones(J-1,1) ], -1:1, J-1, J-1);

    dA = decomposition(A);

    % Inicialización de la matriz de solución U (puntos interiores)
    U_int = zeros(J-1, N+1);
    U_int(:, 1) = u; 

    % Inicialización del error
    errors = zeros(1, N+1);
    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex)) / (max(u_ex) + eps) ;

    for n = 1:N
        u = dA \ (B * u);
        U_int(:, n+1) = u; 
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex)) / (max(u_ex) + eps);
    end
end
