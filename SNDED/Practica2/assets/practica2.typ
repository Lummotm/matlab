#set page(
  paper: "a4",
)
#set math.frac(style: "skewed") //Fuerza que la fracciones sean facheras

= Práctica 2 (Ejercicio 6)
David Nikolov Yordanov

Queremos aproximar la ecuación eliptica:
$ -u'' + u = f(x), " " u(0) = u(1) = 0 $

donde $f(x) = (1 + pi^2) sin(pi x)$.

== La aproximación

Consideramos una función test, $v in V = H^1_0 = { v in S^1: u'(0) = u'(1) = 0 }$.

Multiplicamos pues a la expresión que queremos aproximar:
$ -u''v + u v = f(x) v $

Integrando en $[0,1]$ tenemos:
$ integral_0^1 -u''v + u v " " d x= integral_0^1 f v " " d x $
// Revisar como tagear ecuaciones para poder referenciarlas entre si

Aplicando la regla de la cadena en el primer término de la suma:
$ integral_0^1 -u''v " " d x = underbrace(-u'v]_0^1, 0) + integral_0^1 u'v' " " d x = integral_0^1 u'v " " d x $

Volviendo a la ecuacion primera tenemos que la expresión en forma débil será:
// la primera de la itnegral, no se como referenciarla

$ integral_0^1 u' v' " "d x + integral_0^1 u v " " d x = integral_0^1 f v " " d x $

Consideramos ahora el Método Galerkin, es decir, trabajamos en un espacio discreto $V_h subset V$, como el espacio es finito podemos considerar una base del mismo ${phi_i (x)}_(i=1)^N$.

La solución aproximada se definira como $u_h (x) = sum_(j=1)^N c_j phi_j (x)$. También se puede considerar en vez de la $v$ trabajar con las $phi_i$ de la base.

Llevando la aproximación a la formulación débil nos deja con:
$
  c_j ( integral_0^1 phi_j ' phi_i ' " " d x + integral_0^1 phi_j phi_i " " d x ) = integral_0^1 f phi_i " " d x "    " i = 1,...,N
$

Vectorialmente:
$
  c(K+M) = F
$

Las matrices K y M, en general se denominan, matriz de rigidez y masa respectivamente.

Con coeficientes: $M_(i j) = integral_Omega phi_i phi_j " " d x$, $K_(i j) = integral_Omega phi_i ' phi_j ' " " d x$.

En este trabajo, estudiamos los casos de aproximación por polinomios lineales a trozos y polinomios cuadráticos a trozos.


#pagebreak()

== Elementos Lineleas

Trabajamos en intervalos equiespaciados en el $[0,1]$.

Suponemos el espacio $ V_h = { phi(x) : phi_k |_[x_k,x_(k+1)] "lineal " } $

Consideramos primero un intervalo de referencia $[0,1]$, y definimos las funciones base lineales en este intervalo:
$
  psi_1(x) = 1 - x \
  psi_2(x) = x
$

Que cumplen $psi_1 (0) = psi_2 (1) = 1 ", " psi_1(1) = psi_2(0) = 0$.

Sabemos las considerando el cambio al intervalo de referencia los coeficientes serán de la forma:
$
  M_(i j) = integral_Omega phi_i phi_j " " d x = h integral_0^1 psi_i psi_j " " d x
  K_(i j) = integral_Omega phi_i ' phi_j ' " " d x = h / h^2 integral_0^1 psi_i ' psi_j ' " " d x
$

Donde los $psi_l$ se escogeran de tal forma que mantengan la pendiente, si tiene pendiente negativa se toma $psi_1$, si es positiva entonces se considerara $psi_2$.

Luego tendremos que:
// Revisar si se pueden forzar que las fracciones sean skewed solo para la matriz
$
  M = h/6 mat(
    4, 1, 0, 0, ..., 0, 0;
    1, 4, 1, 0, ..., 0, 0;
    0, 1, 4, 1, ..., 0, 0;
    dots.v, dots.v, dots.v, dots.v, dots.down, dots.v, dots.v;
    0, 0, 0, 0, dots, 1, 4;
  )
$
$
  K = 1/h mat(
    2, -1, 0, 0, ..., 0, 0;
    -1, 2, -1, 0, ..., 0, 0;
    0, -1, 2, -1, ..., 0, 0;
    dots.v, dots.v, dots.v, dots.v, dots.down, dots.v, dots.v;
    0, 0, 0, 0, dots, -1, 2;
  )
$








