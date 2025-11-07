function [Times,Errors,U_int,t,x] = practica1_1(choiceMethod, J_values, N_values, choiceError)
    % Si quiere graficar, dar solo un valor de J, N en general el programa esta diseñado para obtener todas las graficas de una sola vez
    % Dato del problema
    T = 0.5;

    % Condiciones iniciales y solución exacta
    u0_fun = @(x) sin(2*pi*x);
    u_exact = @(x,t) sin(2*pi*x) .* exp(-4*pi^2 * t);

    % Gestión de inputs
    if nargin < 1 || isempty(choiceMethod)
        choiceMethod = 0;
    end

    if nargin < 2 || isempty(J_values)
        J_values = [100];
    end

    if nargin < 3 || isempty(N_values)
        % Condición de estabilidad (mu <= 1/2 asegurado)
        % Uso 2.5 en vez de 2 ya que sino podria pasar qeu mu no este por
        % debajo por errores de redondeo
        N_values = 2.5*J_values.^2*T;
    end

    if nargin < 4 || isempty(choiceError)
        disp("choiceError puede ser:")
        disp("1: Error relativo")
        disp("2: Error máximo")
        disp("Como no se ha introducido nada, usamos error máximo.")
        choiceError = 1; 
    end

    % Exigir el metodo y no crashear directamente si no se pasa nada
    validMethods = 1:3;
    while ~ismember(choiceMethod, validMethods)
        disp("Métodos de resolución:");
        disp("1) Explícito");
        disp("2) Implícito");
        disp("3) Crank-Nicolson");
        choiceMethod = input("Selecciona un método [1-3]: ");
    end

    method_names = {'Explícito', 'Implícito', 'Crank-Nicolson'};


    % Inicialización de los valores de k, h 
    L_j = length(J_values);
    L_n = length(N_values); 
    h_values = 1 ./ J_values;
    k_values = T ./ N_values; 


    % Inicialización de matrices de resultados
    Times = zeros(L_n, L_j);
    Errors = zeros(L_n, L_j);

    % Inicialización de la matriz de la solución (solo se devolverá para un solo J, N)
    U_int = [];

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
            [time, errors, U_int_current] = solver_selection(choiceMethod, params, u0, u_exact, choiceError);

            errmax = max(errors);
            Times(n, j) = time;
            Errors(n, j) = errmax;

            % Guardar la matriz de solución solo si es una única corrida (para visualización)
            if L_n == 1 && L_j == 1
                U_int = U_int_current;
            end
        end
    end

    nombre_metodo = method_names{choiceMethod};

    if  length(N_values) > 1 
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

        % Corregido: removido "\n"
        print("-f1", "k_fijo_var_h_COMBINADO_" + nombre_metodo, "-dpng");

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

        % Corregido: removido "\n"
        print("-f2", "h_fijo_var_k_COMBINADO_" + nombre_metodo, "-dpng");

        % 3) Eficiencia
        figure(3);
        clf;      
        hold on; 

        for n = 1:L_n
            loglog(Times(n,:), Errors(n,:), 'o', 'MarkerSize', 8, 'LineWidth', 2);
        end

        hold off; 
        grid on;
        xlabel('Tiempo de cómputo (s)');
        ylabel('Error máximo');
        title(sprintf('Eficiencia: %s', nombre_metodo));

        legend(legend_text_k, 'Location', 'best');

        % Corregido: removido "\n"
        print("-f3", "eficiencia_" + nombre_metodo, "-dpng")
    end
end



function [time, errors, U_int] = solver_selection(choiceMethod, params, u0, u_exact, choiceError)
    switch choiceMethod
    case 1
        if params.mu <= 1/2 
            tic;
            [U_int, errors] = solve_explicito(params, u0, u_exact, choiceError);
            time = toc;
        else 
            errors(:) = NaN;
            time = NaN;
            U_int = NaN;
        end
    case 2
        tic;
        [U_int, errors] = solve_implicito(params, u0, u_exact, choiceError);
        time = toc;
    case 3
        tic;
        [U_int, errors] = solve_crank(params, u0, u_exact, choiceError);
        time = toc;
    end
end

function [U_int, errors] = solve_explicito(params, u, u_exact, choiceError)
    J = params.J;
    N = params.N;
    mu = params.mu;
    x = params.x;
    t = params.t;

    errors = zeros(1, N+1);

    % Inicialización de la matriz de solución U (interior points)
    U_int = zeros(J-1, N+1);
    U_int(:, 1) = u; % Almacenar la condición inicial

    A = spdiags([mu , (1-2*mu), mu], -1:1, J-1, J-1);

    u_ex = u_exact(x, t(1));
    if choiceError == 1 
        errors(1) = max(abs(u - u_ex)) / (max(u_ex) + eps) ;
    elseif chioceError == 2 
        errors(1) = max(abs(u - u_ex)); 
    end

    for n = 1:N
        u = A * u;
        U_int(:, n+1) = u; % Almacenar el paso de tiempo
        u_ex = u_exact(x, t(n+1));
        if choiceError == 1
            errors(n+1) = max(abs(u - u_ex)) / (max(u_ex) + eps);
        elseif  choiceError == 2
            errors(n+1) = max(abs(u-u_ex));
        end
    end
end

function [U_int, errors] = solve_implicito(params, u, u_exact, choiceError)
    J = params.J;
    N = params.N;
    mu = params.mu;
    x = params.x;
    t = params.t;

    A = spdiags([-mu , (1+2*mu), -mu ], -1:1, J-1, J-1);
    dA = decomposition(A);

    errors = zeros(1, N+1);

    % Inicialización de la matriz de solución U (interior points)
    U_int = zeros(J-1, N+1);
    U_int(:, 1) = u; % Almacenar la condición inicial

    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex)) / max(u_ex);


    for n = 1:N
        u = dA \ u;
        U_int(:, n+1) = u; % Almacenar el paso de tiempo
        u_ex = u_exact(x, t(n+1));
        if choiceError == 1
            errors(n+1) = max(abs(u - u_ex)) / (max(u_ex) + eps);
        elseif  choiceError == 2
            errors(n+1) = max(abs(u-u_ex));
        end
    end
end

function [U_int, errors] = solve_crank(params, u, u_exact, choiceError)
    J = params.J;
    N = params.N;
    mu = params.mu;
    x = params.x;
    t = params.t;

    temp = mu / 2;
    A = spdiags([-temp, (1+mu) , -temp ], -1:1, J-1, J-1);
    B = spdiags([temp,  (1-mu),  temp ], -1:1, J-1, J-1);
    dA = decomposition(A);

    errors = zeros(1, N+1);

    % Inicialización de la matriz de solución U (interior points)
    U_int = zeros(J-1, N+1);
    U_int(:, 1) = u; % Almacenar la condición inicial

    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex)) / max(u_ex);

    for n = 1:N
        u = dA \ (B * u);
        U_int(:, n+1) = u; % Almacenar el paso de tiempo
        u_ex = u_exact(x, t(n+1));
        if choiceError == 1
            errors(n+1) = max(abs(u - u_ex)) / (max(u_ex) + eps);
        elseif  choiceError == 2
            errors(n+1) = max(abs(u-u_ex));
        end
    end
end
