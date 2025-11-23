% Si hice un script para computar los coerfs de cada submatriz M que pasa, no tengo ganas de integrar a mano productos de grado 4
phi_1 = @(x) (2*x-1).*(x-1);
phi_2 = @(x) x.*(1-x)*4;
phi_3 = @(x) (2*x-1).*x;

M_11 = @(x) phi_1(x).*phi_1(x);
M_22 = @(x) phi_2(x).*phi_2(x);
M_33 = @(x) phi_3(x).*phi_3(x);
M_12 = @(x) phi_1(x).*phi_2(x);
M_13 = @(x) phi_1(x).*phi_3(x);
M_23 = @(x) phi_2(x).*phi_3(x);

M(1,1) = integral(M_11,0,1); %= 2/15
M(2,2) = integral(M_22,0,1); %= 8/15
M(3,3) = integral(M_33,0,1); %= 2/15
M(1,2) = integral(M_12,0,1); %= 1/15
M(1,3) = integral(M_13,0,1); %= -1/30
M(2,3) = integral(M_23,0,1); %= 1/15

% Si hice un script para computar los coerfs de cada submatriz K que pasa, no tengo ganas de integrar a mano productos de grado 4
der_phi_1 = @(x) 4*x-3;
der_phi_2 = @(x) 4*(1-2*x);
der_phi_3 = @(x) 4*x-1;

K_11 = @(x) der_phi_1(x).*der_phi_1(x);
K_22 = @(x) der_phi_2(x).*der_phi_2(x);
K_33 = @(x) der_phi_3(x).*der_phi_3(x);
K_12 = @(x) der_phi_1(x).*der_phi_2(x);
K_13 = @(x) der_phi_1(x).*der_phi_3(x);
K_23 = @(x) der_phi_2(x).*der_phi_3(x);

K(1,1) = integral(K_11,0,1); %=7/3
K(2,2) = integral(K_22,0,1); %=16/3
K(3,3) = integral(K_33,0,1); %=7/3
K(1,2) = integral(K_12,0,1); %=-8/3
K(1,3) = integral(K_13,0,1); %=1/3
K(2,3) = integral(K_23,0,1); %=-8/3

disp("Usar (format rational) para verlo en fracciones")
disp("Coefs de M (M es simetrica)")
full(M)
disp("Coefs de K(K es simetrica)")
full(K)

