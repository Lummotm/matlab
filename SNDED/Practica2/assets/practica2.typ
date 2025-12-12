#set page(
  paper: "a4",
)
#let equ(eq, id: none) = {
  let body = if type(id) == none { eq } else if type(id) == label [#eq #id] else [#eq <#id>]
  let numbering = if type(id) != none { "(1)" } else { none }
  set math.equation(numbering: numbering)
  body
}
#set math.equation(supplement: none)

= Práctica 2 (Ejercicio 6)
David Nikolov Yordanov

Queremos aproximar la ecuación eliptica:
$ -u'' + u = f(x), " " u(0) = u(1) = 0 $

donde $f(x) = (1 + pi^2) sin(pi x)$.

== La aproximación

Consideramos una función test, $v in V = H^1_0 = { v in S^1: v(0) = v(1) = 0 }$.

Multiplicamos pues a la expresión que queremos aproximar:
$ -u''v + u v = f(x) v $

Integrando en $[0,1]$ tenemos:
#equ($ integral_0^1 -u''v + u v " " d x= integral_0^1 f v " " d x $, id: <eq:eq1>)
// Revisar como tagear ecuaciones para poder referenciarlas entre si

Aplicando la regla de la cadena en el primer término de la suma:
$ integral_0^1 -u''v " " d x = underbrace([-u'v]_0^1, 0) + integral_0^1 u'v' " " d x = integral_0^1 u'v " " d x $

Volviendo a la ecuación @eq:eq1 tenemos que la expresión en forma débil será:

$ integral_0^1 u' v' " "d x + integral_0^1 u v " " d x = integral_0^1 f v " " d x $

Consideramos ahora el Método Galerkin, es decir, trabajamos en un espacio discreto $V_h subset V$, como el espacio es finito podemos considerar una base del mismo ${phi_i (x)}_(i=1)^N$.

La solución aproximada se definira como $u_h (x) = sum_(j=1)^N c_j phi_j (x)$. Como la funciónn base se puede expresar como suma de coeficientes de las $phi_i$ entonces trabajamos con estas en vez de $v$.

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

== Elementos Lineales

Trabajamos en intervalos equiespaciados en el $[0,1]$.

Suponemos el espacio $ V_h = { phi(x) : phi_k |_[x_k,x_(k+1)] "lineal " } $

Para cada elemento $[x_(k-1),x_(k+1)]$ podemos definir los polinomios como:
$
  phi_j (x) := cases(
    (x - x_(j-1)) / (x_(j) - x_(j-1)) "si" x_(j-1) <= x <= x_(j),
    (x_(j+1) - x) / (x_(j+1) - x_j) "si" x_(j) <= x <= x_(j+1),
  )
$

Consideramos primero un intervalo de referencia $[0,1]$, y definimos las funciones base lineales en este intervalo:
$
      psi_1(zeta) = zeta "  " & psi_1 ' (zeta) = 1 \
  psi_2(zeta) = 1 - zeta "  " & psi_2 ' (zeta) = -1
$

Que cumplen $psi_1 (0) = psi_2 (1) = 0 ", " psi_1(1) = psi_2(0) = 1$.

El cambio de parametro será: $x(zeta) = x_j + zeta h$, luego $d x = h d zeta$.

=== Matrices de masa y rigidez

Sabemos las considerando el cambio al intervalo de referencia los coeficientes serán de la forma:
// falta pasar bien la sección del calculo esta en la tablet
$
  M_(i j) = integral_Omega phi_i phi_j " " d x = integral_0^1 psi_i psi_j " " h d zeta =h integral_0^1 psi_i psi_j " " d zeta
$
$
  K_(i j) = integral_Omega phi_i ' phi_j ' " " d x = integral_0^1 psi_i ' 1/h " " psi_j ' " " 1/h " " h " " d zeta = 1/h integral_0^1 psi_i ' psi_j ' d zeta
$

Consideramos pues las diferentes posiciones principales para $M$, es decir:
$
  M_(j j) & = integral_(x_(j-1))^(x_(j+1)) phi_j^2 d x
            = integral_(x_(j-1))^x_j phi_j^2 d x + integral_(x_(j))^x_(x_(j+1)) phi_j^2 d x \
          & = underbrace(integral_0^1 psi_1^2 h " " d zeta, "parte que baja") + underbrace(integral_0^1 psi_2^2 h " " d zeta, "parte que sube")
            = underbrace(2 h integral_0^1 (1 - zeta)^2 d zeta, "simetría de las funciones") =(2h)/3
$
$
  M_(j j-1) = integral_x_(j-1)^x_j phi_j phi_(j-1) d x
  = integral_0^1 underbrace((1-zeta), phi_j "baja") underbrace(zeta, phi_(j-1) "sube") h " " d zeta
  = h/6
$
$
  M_(j j+1) = integral_x_(j-1)^x_j phi_j phi_(j+1) d x
  = integral_0^1 underbrace(zeta, phi_(j) "sube")
  underbrace((1-zeta), phi_(j+1) "baja") h " " d zeta
  = h/6
$

#pagebreak()

Las derivadas de las funciones $psi_i$ son o $1$ o $-1$.

Consideramos ahora las posiciones de $K$, es decir:

$
  K_(j j) & = integral_(x_(j-1))^(x_(j+1)) (phi_j ')^2 d x
            = underbrace(integral_0^1 (1/h)^2 " " h " " d zeta, "parte que sube")
            + underbrace(integral_0^1 (-1/h)^2 h " "d zeta, "parte que baja")
            = 2 1/h integral_0^1 d zeta = 2/h
$
$
  K_(j j-1) & = underbrace(integral_(x_(j-1))^(x_(j)) phi_(j-1) ' phi_j ' d x, "no nulo en" [x_(j-1),x_j]) = integral_0^1 -1/h dot 1/h " " h d zeta = -1/h = -1/h
$

$
  K_(j j+1) & = underbrace(integral_(x_(j+1))^(x_j) phi_(j) ' phi_(j+1) ' d x, "no nulo un " [x_j, x_(j+1)]) = integral_0^1 -1/h dot 1/h " " h d zeta = -1/h
$

Entonces las matrices serán de la forma:
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

=== Cálculo del Vector de Carga (F)

El término de la fuerza $F$ se calcula ensamblando las contribuciones locales.

$
  F_i = integral_0^1 f(x(zeta)) phi_i (x(zeta)) h d zeta
$


#pagebreak()

== Elementos Cuadráticos

Para mejorar el orden de convergencia, consideramos el espacio de elementos finitos cuadráticos a trozos. En el elemento de referencia $[0,1]$, utilizamos tres funciones de forma (Lagrange de orden 2):

$
  psi_1(zeta) = (1-zeta)(1-2zeta) "   " psi_1 '(zeta) & = 4 zeta - 3 \
      psi_2(zeta) = 4zeta(1-zeta) "   " psi_2 '(zeta) & = 4 ( 1 - 2 zeta) \
       psi_3(zeta)= zeta(2zeta-1) "   " psi_3 '(zeta) & = 4 zeta - 1
$

Estas funciones satisfacen $psi_1(0)=1$, $psi_2(1/2)=1$ y $psi_3(1)=1$.

=== Matrices Elementales Cuadráticas

Al integrar estas funciones en $[0,1]$, obtenemos matrices locales de dimensión $3 times 3$.

*Matriz de Masa Local ($M^e$):*
$
  M^e = h/30 mat(
    4, 2, -1;
    2, 16, 2;
    -1, 2, 4
  )
$

*Matriz de Rigidez Local ($K^e$):*
$
  K^e = 1/(3h) mat(
    7, -8, 1;
    -8, 16, -8;
    1, -8, 7
  )
$

El ensamblaje produce matrices globales se realiza sumando las contribuciones de las matrices en los respectivos nodos compartidos por elemento.

#pagebreak()
== Resultados Numéricos

Se ha implementado el método en MATLAB para mallas uniformes. Un aspecto crítico observado durante la experimentación es la influencia de la integración numérica en el cálculo del vector de carga $F$.

=== Comparación de Cuadraturas

Se probaron tres métodos de integración para calcular $integral f phi_i$: Rectángulo, Punto Medio y Simpson. A continuación se muestran las pendientes de convergencia ($p$) obtenidas experimentalmente ($E approx C h^p$) para cada caso:

#figure(
  table(
    columns: (auto, auto, auto),
    inset: 10pt,
    align: center,
    [*Método de Integración*], [*Pendiente Lineal ($p$)*], [*Pendiente Cuadrática ($p$)*],
    [Rectángulo (Izq)], [2.00], [2.00],
    [Punto Medio], [2.00], [2.01],
    [Simpson], [2.00], [3.00],
  ),
  caption: [Orden de convergencia observado según la regla de cuadratura.],
)

*Análisis:*
- Para *elementos lineales*, el error teórico es $O(h^2)$. Cualquier cuadratura de orden al menos 1 o 2 es suficiente, por lo que todos los métodos arrojan la pendiente correcta $p approx 2$.
- Para *elementos cuadráticos*, el error teórico es $O(h^3)$. Sin embargo, al usar reglas de bajo orden (Rectángulo o Punto Medio), el error de integración numérica es $O(h^2)$. Este error domina sobre el error de aproximación, degradando la convergencia global a $p approx 2$.
- Solo al utilizar la *Regla de Simpson* (orden superior), el error de integración se vuelve despreciable frente al de aproximación, permitiendo recuperar el orden óptimo $p approx 3$.

=== Tabla de Convergencia (Usando Simpson)

A continuación se detallan los errores absolutos utilizando la cuadratura adecuada (Simpson).

#figure(
  table(
    columns: (auto, auto, auto, auto),
    inset: 10pt,
    align: center,
    [*N*], [*h*], [*Error Lineal* ($L^2$)], [*Error Cuadrático* ($L^2$)],
    [10], [0.1], [5.89e-3], [1.27e-4],
    [20], [0.05], [1.47e-3], [1.58e-5],
    [40], [0.025], [3.68e-4], [1.97e-6],
    [80], [0.0125], [9.20e-5], [2.46e-7],
  ),
  caption: [Evolución del error $L^2$.],
)

El análisis logarítmico con estos datos confirma que:
1. Los elementos lineales convergen como $O(h^2)$.
2. Los elementos cuadráticos convergen como $O(h^3)$.
