function [Times,Errors] = practica1_1(choiceMethod, h_values, k_values)
    if nargin < 1 || isempty(choiceMethod)
        choiceMethod = 0;
    end

    % Si no se introduce nada manejo los inputs por mi cuenta
    if nargin < 2 || isempty(h_values)
        h_values = 10.^(-(1:4)); 
    end

    if nargin < 3 || isempty(k_values)
        k_values = h_values; % Por defecto, k = h
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

    % Dato del problema
    T = 0.5;

    % Condiciones iniciales y solución exacta
    u0_fun = @(x) sin(pi*x);
    u_exact = @(x,t) sin(pi*x) .* exp(-pi^2 * t);

    % Inicialización de matrices de resultados
    L = length(k_values);
    Times = zeros(L, L);
    Errors = zeros(L, L);

    for i = 1:L
        k = k_values(i);

        for m = 1:L
            h = h_values(m);
            mu = k / h^2;
            J = ceil(1 / h);
            N = ceil(T / k);

            % Malla
            t = linspace(0, T, N+1)';
            x = linspace(0, 1, J+1)';
            x = x(2:end-1); 

            % Condición inicial
            u0 = u0_fun(x);

            % Resolver
            [time, errors, ~] = solver_selection(choiceMethod, J, N, mu, u0, x, t, u_exact);

            errmax = max(errors);
            Times(i, m) = time;
            Errors(i, m) = errmax;
        end
    end

    nombre_metodo = method_names{choiceMethod};

    % 1) k fijo, variando h
    for i = 1:L
        subplot(ceil(L/2), 2, i);
        loglog(Times(i,:), Errors(i,:), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
        grid on;
        xlabel('Tiempo de ejecución');
        ylabel('Error máximo');
        title(sprintf('k = %.0e fijo, variando h', k_values(i)));
    end

    % 2) h fijo, variando k
    for m = 1:L
        subplot(ceil(L/2), 2, m);
        loglog(Times(:,m), Errors(:,m), 's-', 'LineWidth', 2, 'MarkerSize', 8);
        grid on;
        xlabel('Tiempo de ejecución');
        ylabel('Error máximo');
        title(sprintf('h = %.0e fijo, variando k', h_values(m)));
    end

    % 3) Eficiencia
    loglog(Times(:), Errors(:), 'o', 'MarkerSize', 8);
    grid on;
    xlabel('Tiempo de ejecución');
    ylabel('Error máximo');
    title(sprintf('Eficiencia: %s', nombre_metodo));

end

function [time, errors, u] = solver_selection(choice, j, n, mu, u0, x, t, u_exact)
    switch choice
    case 1
        if mu <= 1/2 
            tic;
            [u, errors] = solve_explicito(j, n, mu, u0, x, t, u_exact);
            time = toc;
        else 
            errors(:) = NaN;
            time = NaN;
            u(:) = NaN;
        end
    case 2
        tic;
        [u, errors] = solve_implicito(j, n, mu, u0, x, t, u_exact);
        time = toc;
    case 3
        tic;
        [u, errors] = solve_crank(j, n, mu, u0, x, t, u_exact);
        time = toc;
    end
end

function [u, errors] = solve_explicito(j, n, mu, u, x, t, u_exact)
    errors = zeros(1, n+1);

    A = spdiags([mu , (1-2*mu), mu], -1:1, j-1, j-1);

    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex));

    for n = 1:n
        u = A * u;
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex));
    end
end

function [u, errors] = solve_implicito(j, n, mu, u, x, t, u_exact)
    A = spdiags([-mu , (1+2*mu) , -mu ], -1:1, j-1, j-1);
    dA = decomposition(A);

    errors = zeros(1, n+1);

    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex));

    for n = 1:n
        u = dA \ u;
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex));
    end
end

function [u, errors] = solve_crank(j, n, mu, u, x, t, u_exact)
    temp = mu / 2;
    A = spdiags([-temp, (1+mu) , -temp ], -1:1, j-1, j-1);
    B = spdiags([temp  (1-mu)  temp ], -1:1, j-1, j-1);
    dA = decomposition(A);

    errors = zeros(1, n+1);

    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex));

    for n = 1:n
        u = dA \ (B * u);
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex));
    end
end
