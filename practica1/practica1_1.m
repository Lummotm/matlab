function practica1_1(choice,J)

    % Función auxiliar para verificar que un valor es entero
    % Fuente: https://www.analyticslane.com/2022/02/09/comprobar-si-un-valor-es-entero-en-matlab/
    isint = @(x) round(x) == x;

    % Comprobación de inputs opcionales
    if nargin < 1 || isempty(choice)
        choice = 0;
    end
    if nargin < 2 || isempty(J)
        J = input("Introduce el número de nodos: ");
        % Verificar que J sea un entero positivo
        while ~isint(J) || J <= 0
            disp("El valor de J debe ser un entero positivo.")
            J = input("Introduce el número de nodos: ");
        end
    end

    % Opciones válidas para el método
    validOptions = [1:3];
    % Pedir método hasta que se introduzca una opción válida
    while ~ismember(choice,validOptions)
        disp("=== Métodos de resolución ===");
        disp("1) Explícito");
        disp("2) Implícito");
        disp("3) Crank-Nicolson");
        choice = input("Selecciona un método [1-3]: ");
    end


    % Falta hacer la casuistica, teniendo en cuenta valores de k fijos y variando h y a la inversa
    % Datos del problema
    T = 0.5;
    % Definición de parámetros de la malla para estabilidad
    h = 1/J;             % paso espacial
    N = ceil(2*T/h^2);   % número de pasos de tiempo para cumplir mu <= 1/2
    k = T/N;             % paso temporal
    % Definición de la malla espacial
    x = linspace(0,1,J+1)';
    x = x(2:end-1);      % nodos internos (excluyendo contorno)
    t = linspace(0,T,N+1)';
    % Condición inicial (autofunción del problema)
    u0_fun = @(x) sin(pi*x);
    u0 = u0_fun(x);
    % Parámetro de estabilidad
    mu = k/h^2;

    % x sera los nodos de x, mientras que t sera un t para cada tiempo
    u_exact = @(x,t) (sin(pi*x) .* exp(-pi^2 *t));

    % Llamada al solver correspondiente
    switch choice
    case 1
        tic
        [u, errors] = solve_explicito(J,N,mu,u0,x,t,u_exact);
        toc
    case 2
        tic
        [u, errors] = solve_implicito(J,N,mu,u0,x,t,u_exact);
        toc
    case 3
        tic
        [u, errors] = solve_crank(J,N,mu,u0,x,t,u_exact);
        toc
    end

    fprintf('Error en t=0: %.6e\n', errors(1));
    fprintf('Error en t=T: %.6e\n', errors(end));
    fprintf('Error máximo: %.6e\n', max(errors(2:end)));
    fprintf('Error mínimo: %.6e\n', min(errors(2:end)));

    figure;
    loglog(t, errors, 'o-');
    xlabel('Tiempo t');
    ylabel('Error máximo');
    title('Error máximo vs tiempo (log-log)');
    grid on;
end

function [u, errors] = solve_explicito(J,N,mu,u,x,t,u_exact)
    % Defino matriz dispersa del problema
    A = spdiags([mu 1-2*mu mu],-1:1,J-1,J-1);

    % Vector para almacenar errores
    errors = zeros(1,N+1);

    % Error en t=0
    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex));

    % Iteramos en el cada orden de tiempo
    % En el metodo explicito las iteraciones son de la forma u_n+1 = A * u_n
    for n = 1:N
        u = A * u;
        % Calcular error en este paso de tiempo
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex));
    end

    u = [0; u; 0];
end

function [u, errors] = solve_implicito(J,N,mu,u,x,t,u_exact)
    % Defino matriz dispersa del problema
    A = spdiags([-mu 1+2*mu -mu],-1:1,J-1,J-1);
    dA = decomposition(A);

    % Vector para almacenar errores
    errors = zeros(1,N+1);

    % Error en t=0
    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex));

    for n = 1:N
        u = dA \ u;
        % Calcular error en este paso de tiempo
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex));
    end

    u = [0; u; 0];
end

function [u, errors] = solve_crank(J,N,mu,u,x,t,u_exact)
    temp = mu/2;
    A = spdiags([-temp 1+mu -temp],-1:1,J-1,J-1);
    B = spdiags([-temp 1-mu -temp],-1:1,J-1,J-1);   
    dA = decomposition(A);

    % Vector para almacenar errores
    errors = zeros(1,N+1);

    % Error en t=0
    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex));

    for n = 1:N
        u = dA \ (B * u);
        % Calcular error en este paso de tiempo
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex));
    end

    u = [0; u; 0];
end
