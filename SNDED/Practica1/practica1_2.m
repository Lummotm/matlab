function [Times, Errors, U_int, t, x] = practica1_2(choiceMethod, J_values, N_values, T)
    % Resuelve la ecuación del calor 1D NO HOMOGÉNEA (u_t = u_xx + f(x,t))
    % con condiciones de borde no homogéneas.
    %
    % Input: (choiceMethod, J_values, N_values)
    % - choiceMethod: metodo escogido 
    %       - 1) Explícito
    %       - 2) Implícito
    %       - 3) Crank-Nicolson
    % - J_values: valores de J. Vector o escalar.
    % - N_values: valores de N. Vector o escalar. Si se omite ([]),
    %             se calcula N para mu=0.4 (régimen estable).
    % - T: valor final de tiempo.
    %
    % Output:
    % - Times: Matriz (L_n x L_j) de tiempos de cómputo.
    % - Errors: Matriz (L_n x L_j) de errores máximos relativos.
    % - U_int, t, x: Solución y mallas (solo si J y N son escalares).

    if nargin < 1 || isempty(choiceMethod)
        choiceMethod = 0; % Valor por defecto si no se introduce
    end
    
    if nargin < 2 || isempty(J_values)
        J_values = [100]; 
    end

     % En la practica1_2 debería de variar el valor de T
    if isempty(T)
        T = 5;
    end

    if nargin < 3 || isempty(N_values)
        % Si no se da, se calcula N para CUMPLIR LA CONDICIÓN DE ESTABILIDAD
        % del método explícito (mu <= 1/2).
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

   

    % Condiciones iniciales y solución exacta 
    u_exact = @(x,t) t .* cos(pi*x/2) .^ 2 ./ (t+1);

    % u0 sera la dada por el exacto en el borde t=0
    u0_fun = @(x) u_exact(x,0);
    f_fun = @(x,t) (cos(pi*x/2) ./ (t+1)) .^ 2 + (pi^2*t .* cos(pi*x)) ./ (2*(t+1));
    
    % Condición borde (0,t)
    b_fun = @(t) t ./ (t+1);


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
            [time, errors, U_int_current] = solver_selection(choiceMethod, params, u0, u_exact, f_fun, b_fun);

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
        print("-f1", "P2_k_fijo_var_h_" + nombre_metodo, "-dpng");

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
        print("-f2", "P2_h_fijo_var_k_" + nombre_metodo, "-dpng");
    end
end


function [time, errors, U_int] = solver_selection(choiceMethod, params, u0, u_exact, f_fun, b_fun)
    J = params.J;
    N = params.N;
    
    switch choiceMethod
    case 1 % Explícito
        % Requiere condición: mu <= 1/2 para estabilidad
        if params.mu <= 1/2 
            tic;
            [U_int, errors] = solve_explicito(params, u0, u_exact, f_fun, b_fun); 
            time = toc;
        else 
            % Inestable: devolvemos NaN
            errors = NaN(1, N + 1);
            time = NaN;
            U_int = NaN(J - 1, N + 1);
        end
    case 2 % Implícito 
        tic;
        [U_int, errors] = solve_implicito(params, u0, u_exact, f_fun, b_fun);
        time = toc;
    case 3 % Crank-Nicolson 
        tic;
        [U_int, errors] = solve_crank(params, u0, u_exact, f_fun, b_fun);
        time = toc;
    end
end

function [U_int, errors] = solve_explicito(params, u, u_exact, f_fun, b_fun)
    J = params.J; 
    N = params.N; 
    mu = params.mu;
    x = params.x; 
    t = params.t; 
    k = t(2) - t(1);

    A = spdiags([mu*ones(J-1,1) , (1-2*mu)*ones(J-1,1), mu*ones(J-1,1)], -1:1, J-1, J-1);

    U_int = zeros(J-1, N+1); 
    U_int(:, 1) = u; 
    errors = zeros(1, N+1);
    
    % Error inicial (t=0)
    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex)) / (max(u_ex) + eps);
    b = zeros(J-1 ,1);

    for n = 1:N
        % Evaluación fuente y borde en t_n 
        f = f_fun(x,t(n));  
        b(1) = mu * b_fun(t(n));

        % Cálculo
        u = A * u + k * f + b;

        % Guardamos los datos de la iteración
        U_int(:, n+1) = u; 
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex)) / (max(u_ex) + eps);
    end
end

function [U_int, errors] = solve_implicito(params, u, u_exact, f_fun, b_fun)
    J = params.J; 
    N = params.N; 
    mu = params.mu;
    x = params.x; 
    t = params.t; 
    k = t(2) - t(1);

    A = spdiags([-mu*ones(J-1,1) , (1+2*mu)*ones(J-1,1), -mu*ones(J-1,1) ], -1:1, J-1, J-1);
    
    % Pre-factorización para optimizar resolución en bucle
    dA = decomposition(A);

    U_int = zeros(J-1, N+1);
    U_int(:, 1) = u; 
    errors = zeros(1, N+1);
    
    % Error incial
    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex)) / (max(u_ex) + eps);
    b = zeros(J-1, 1);

    for n = 1:N
        % Evaluación fuente y borde en t_{n+1}
        f = f_fun(x, t(n+1));
        b(1) = mu * b_fun(t(n+1));
        
        % Cálculo
        LD = u + k * f + b; % Lado derecho
        u = dA \ LD;

        % Guardamos los datos de la iteración
        U_int(:, n+1) = u; 
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex)) / (max(u_ex) + eps);
    end
end

function [U_int, errors] = solve_crank(params, u, u_exact, f_fun, b_fun)
    J = params.J; 
    N = params.N; 
    mu = params.mu;
    x = params.x; 
    t = params.t; 
    k = t(2) - t(1);

    temp = mu / 2;
    A = spdiags([-temp*ones(J-1,1), (1+mu)*ones(J-1,1) , -temp*ones(J-1,1) ], -1:1, J-1, J-1);
    B = spdiags([temp*ones(J-1,1),  (1-mu)*ones(J-1,1),  temp*ones(J-1,1) ], -1:1, J-1, J-1);
    
    dA = decomposition(A);

    U_int = zeros(J-1, N+1); 
    U_int(:, 1) = u; 
    errors = zeros(1, N+1);
    
    % Error inicial
    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex)) / (max(u_ex) + eps);
    
    % Asignamos primero los valores en t_n
    f_prev = f_fun(x, t(1)); 
    b_prev = zeros(J-1, 1); 
    b_curr = zeros(J-1, 1);
    b_prev(1) = temp * b_fun(t(1)); 

    for n = 1:N
        % Evaluación fuente y borde en t_{n+1}
        f_curr = f_fun(x, t(n+1));
        b_curr(1) = temp * b_fun(t(n+1)); 

        LD = (B*u) + (k/2) * (f_prev + f_curr) + b_prev + b_curr; % Lado derecho
        u = dA \ LD;

        % Actualizar estados previos
        b_prev(1) = b_curr(1); 
        f_prev = f_curr;

        % Guardamos los datos de la iteración
        U_int(:, n+1) = u; 
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex)) / (max(u_ex) + eps);
    end
end
