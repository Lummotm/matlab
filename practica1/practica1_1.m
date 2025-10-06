function practica1_1(choiceMethod)
    % Error máximo (se supone que es un problema de precisión de matlab)
    EMAX=1e18;

    % Comprobación de inputs opcionales
    if nargin < 1 || isempty(choiceMethod)
        choiceMethod = 0;
    end

    validMethods = 1:3;
    while ~ismember(choiceMethod, validMethods)
        disp("=== Métodos de resolución ===");
        disp("1) Explícito");
        disp("2) Implícito");
        disp("3) Crank-Nicolson");
        choiceMethod = input("Selecciona un método [1-3]: ");
    end

    % Dato del problema
    T = 0.5;

    % Condiciones iniciales 
    u0_fun = @(x) (sin(pi*x)); % autofunción 

    % Función de solución 
    u_exact = @(x,t) (sin(pi*x) .* exp(-pi^2 *t));


    % Indice de iteracion (valores de h y k que van a haber)
    L = 5;


    % Hago print de la primera fila en la tabla
    fprintf("  k     |   h   ")
    %fprintf("k|h\t")
    for l = 1:L
        H(l) = 10^(-l);
        fprintf("|%.6e\t",H(l));
    end
    K = H;

    for i = 1:L 
        k = K(i);
        fprintf("\n%.6e\t",k);
        for m = 1:L 
            h = H(m);  % Variamos h
            mu = k / h^2;
            J = ceil(1 / h);       % nodos totales
            N = ceil(T / k);       % pasos de tiempo

            % Definición de la malla
            t = linspace(0,T,N+1)';
            x = linspace(0,1,J+1)';
            x = x(2:end-1);        % nodos internos

            % Evaluación de la condición inicial
            u0 = u0_fun(x);

            % Usamos tantas variables dentro de la función porque no repetir cuentas en cada uno de los métodos
            [time, errors] = solver_selection(choiceMethod,J,N,mu,u0,x,t,u_exact);

            errmax = max(errors);
            if isinf(errmax) || isnan(errmax) || errmax >= EMAX
                % Si el error se pasa de un error logico escribimos *** en vez del error
                fprintf("|************\t")            
            else
                fprintf("|%.6e\t",errmax);           
            end
            times(:,i) = time;
        end    
    end
    fprintf("\n")
end

function [time,errors] = solver_selection(choice,j,n,mu,u0,x,t,u_exact)
    % llamada al solver correspondiente
    switch choice
    case 1
        tic
        [~, errors] = solve_explicito(j,n,mu,u0,x,t,u_exact);
        time = toc;
    case 2
        tic
        [~, errors] = solve_implicito(j,n,mu,u0,x,t,u_exact);
        time = toc;
    case 3
        tic
        [~, errors] = solve_crank(j,n,mu,u0,x,t,u_exact);
        time = toc;
    end
end

function [u, errors] = solve_explicito(j,n,mu,u,x,t,u_exact)
    errors = zeros(1,n+1);    
    if mu >= 1/2    
        errors(1) = 1e30;
    else
        % defino matriz dispersa del problema
        A = spdiags([mu 1-2*mu mu],-1:1,j-1,j-1);

        % vector para almacenar errores
        errors = zeros(1,n+1);

        % error en t=0
        u_ex = u_exact(x, t(1));
        errors(1) = max(abs(u - u_ex));

        % iteramos en cada paso de tiempo
        for n = 1:n
            u = A * u;
            % calcular error en este paso de tiempo
            u_ex = u_exact(x, t(n+1));
            errors(n+1) = max(abs(u - u_ex));
        end
    end

end

function [u, errors] = solve_implicito(j,n,mu,u,x,t,u_exact)
    % defino matriz dispersa del problema
    A = spdiags([-mu 1+2*mu -mu],-1:1,j-1,j-1);
    dA = decomposition(A);

    % vector para almacenar errores
    errors = zeros(1,n+1);

    % error en t=0
    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex));

    for n = 1:n
        u = dA \ u;
        % calcular error en este paso de tiempo
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex));
    end
end

function [u, errors] = solve_crank(j,n,mu,u,x,t,u_exact)
    temp = mu/2;
    A = spdiags([-temp 1+mu -temp],-1:1,j-1,j-1);
    B = spdiags([-temp 1-mu -temp],-1:1,j-1,j-1);   
    dA = decomposition(A);

    % vector para almacenar errores
    errors = zeros(1,n+1);

    % error en t=0
    u_ex = u_exact(x, t(1));
    errors(1) = max(abs(u - u_ex));

    for n = 1:n
        u = dA \ (B * u);
        % calcular error en este paso de tiempo
        u_ex = u_exact(x, t(n+1));
        errors(n+1) = max(abs(u - u_ex));
    end
end
