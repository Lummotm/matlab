% Esto equivale a
phi =  {
    @(x) (2*x-1).*(x-1);
    @(x) x.*(1-x)*4;
    @(x) (2*x-1).*x;
    };
% esto:
% phi{1} = @(x) (2*x-1).*(x-1);
% phi{2} = @(x) x.*(1-x)*4;
% phi{3} = @(x) (2*x-1).*x;

der_phi = {
    @(x) 4*x-3;
    @(x) 4*(1-2*x);
    @(x) 4*x-1;
    };

N = length(phi);

M = zeros(N,N);
K = zeros(N,N);

for i = 1:N
    for j = 1:N
        producto = @(x) phi{i}(x) .* phi{j}(x);
        M(i,j) = integral(producto,0,1);
        producto_derivada = @(x) der_phi{i}(x) .* der_phi{j}(x);
        K(i,j) = integral(producto_derivada,0,1);
    end
end

format rational % Forzar el uso de formato racional para que se vea feten

disp("Coefs de M (M es simetrica)")
disp(M)
disp("Coefs de K(K es simetrica)")
disp(K)

