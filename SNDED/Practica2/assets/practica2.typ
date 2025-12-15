#set page(
  paper: "a4",
)

#set text(
  lang: "es",
  font: "New Computer Modern Math",
)

// Solo numerar ecuaciones especificas (sacado de la documentación)

#let equ(eq, id: none) = {
  let body = if type(id) == none { eq } else if type(id) == label [#eq #id] else [#eq <#id>]
  let numbering = if type(id) != none { "(1)" } else { none }
  set math.equation(numbering: numbering)
  body
}

#set math.equation(supplement: none)

// Portada
#align(center, [
  #v(5cm)

  // Título del Trabajo
  #box(text(size: 20pt, weight: "bold", [Práctica 2 (Elementos Finitos)]))
  #v(0.2cm)
  #box(text(size: 14pt, [David Nikolov Yordanov]))
  #v(0.2cm)
  #box(text(size: 12pt, [15/12/2025]))

  #v(4cm)
  // Índice
  #set align(left)

  #outline(
    title: "Índice",
    depth: 4,
  )
])

#pagebreak()

== Introducción

Queremos aproximar la ecuación eliptica:
$ -u'' + u = f(x), " " u(0) = u(1) = 0 $

donde $f(x) = (1 + pi^2) sin(pi x)$.

La solución analítica exacta de esta ecuación diferencial, con las condiciones de contorno dadas y la función fuente $f(x)$, es:
$ u(x) = sin(pi x) $

== Discretizado del espacio

Consideramos una función test, $v in V = H^1_0 = { v in S^1: v(0) = v(1) = 0 }$.

Multiplicamos pues a la expresión que queremos aproximar:
$ -u''v + u v = f(x) v $

Integrando en $[0,1]$ tenemos:
#equ($ integral_0^1 -u''v + u v " " d x= integral_0^1 f v " " d x $, id: <eq:eq1>)
// Revisar como tagear ecuaciones para poder referenciarlas entre si

Aplicando la regla de la cadena en el primer término de la suma:
$ integral_0^1 -u''v " " d x = underbrace([-u'v]_0^1, 0) + integral_0^1 u'v' " " d x = integral_0^1 u'v " " d x $

Volviendo a la ecuación *(@eq:eq1)* tenemos que la expresión en forma débil será:

$ integral_0^1 u' v' " "d x + integral_0^1 u v " " d x = integral_0^1 f v " " d x $

Consideramos ahora el Método Galerkin, es decir, trabajamos en un espacio discreto $V_h subset V$, como el espacio es finito podemos considerar una base del mismo ${phi_i (x)}_(i=1)^N$.

La solución aproximada se definirá como $u_h (x) = sum_(j=1)^N c_j phi_j (x)$. Como la función base se puede expresar como suma de coeficientes de las $phi_i$ entonces trabajamos con estas en vez de $v$.

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
      psi_1(z) = z "  " & psi_1 ' (z) = 1 \
  psi_2(z) = 1 - z "  " & psi_2 ' (z) = -1
$

Que cumplen $psi_1 (0) = psi_2 (1) = 0 ", " psi_1(1) = psi_2(0) = 1$.

El cambio de parametro será: $x(z) = x_j + z h$, luego $d x = h d z$.

=== Matrices de Masa ($M$) y Rigidez ($K$)

Considerando el cambio al intervalo de referencia, los coeficientes serán:
// falta pasar bien la sección del calculo esta en la tablet
$
  M_(i j) = integral_Omega phi_i phi_j " " d x = integral_0^1 psi_i psi_j " " h d z =h integral_0^1 psi_i psi_j " " d z
$
$
  K_(i j) = integral_Omega phi_i ' phi_j ' " " d x = integral_0^1 psi_i ' 1/h " " psi_j ' " " 1/h " " h " " d z = 1/h integral_0^1 psi_i ' psi_j ' d z
$

Consideramos pues las diferentes posiciones principales para $M$, es decir:
$
  M_(j j) & = integral_(x_(j-1))^(x_(j+1)) phi_j^2 d x
            = integral_(x_(j-1))^x_j phi_j^2 d x + integral_(x_(j))^(x_(j+1)) phi_j^2 d x \
          & = underbrace(integral_0^1 psi_1^2 h " " d z, "parte que baja") + underbrace(integral_0^1 psi_2^2 h " " d z, "parte que sube")
            = underbrace(2 h integral_0^1 (1 - z)^2 d z, "simetría de las funciones") =(2h)/3
$
$
  M_(j j-1) = integral_x_(j-1)^x_j phi_j phi_(j-1) d x
  = integral_0^1 underbrace((1-z), phi_j "baja") underbrace(z, phi_(j-1) "sube") h " " d z
  = h/6
$
$
  M_(j j+1) = integral_x_(j-1)^x_j phi_j phi_(j+1) d x
  = integral_0^1 underbrace(z, phi_(j) "sube") " "
  underbrace((1-z), phi_(j+1) "baja") h " " d z
  = h/6
$

#pagebreak()

Las derivadas de las funciones $psi_i$ son $1$ o $-1$.

Consideramos ahora las posiciones de $K$, es decir:

$
  K_(j j) & = integral_(x_(j-1))^(x_(j+1)) (phi_j ')^2 d x
            = underbrace(integral_0^1 (1/h)^2 " " h " " d z, "parte que sube")
            + underbrace(integral_0^1 (-1/h)^2 h " "d z, "parte que baja")
            = 2 1/h integral_0^1 d z = 2/h
$
$
  K_(j j-1) & = underbrace(integral_(x_(j-1))^(x_(j)) phi_(j-1) ' phi_j ' d x, "no nulo en" [x_(j-1),x_j]) = integral_0^1 -1/h dot 1/h " " h d z = -1/h
$

$
  K_(j j+1) & = underbrace(integral_(x_(j+1))^(x_j) phi_(j) ' phi_(j+1) ' d x, "no nulo un " [x_j, x_(j+1)]) = integral_0^1 -1/h dot 1/h " " h d z = -1/h
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

=== Cálculo del Vector de Carga ($F$)

El término de la fuerza $F$ se calcula ensamblando las contribuciones locales.

$
  F_i = integral_0^1 f(x(z)) phi_i (x(z)) h d z
$


#pagebreak()

== Elementos Cuadráticos

Trabajamos ahora con los polinomios cuadráticos a trozos.

A diferencia de los elementos lineales, introducimos un nodo interno en cada elemento. Esto implica que las funciones de base asociadas a los vértices tienen soporte en 2 elementos (como en el caso lineal), mientras que las de los nodos centrales tienen soporte únicamente en el propio elemento.

El tamaño del elemento es $h = x_(2k+1) - x_(2k-1)$, que es uniforme para todos los elementos al ser la malla equiespaciada.

Utilizamos el intervalo de referencia $[0,1]$ y las funciones de forma cuadráticas (polinomios de Lagrange):

$
  phi_1(z) = (2z-1)(z-1) = 2z^2 - 3z + 1 \
  phi_2(z) = 4z(1-z) = 4z - 4z^2 \
  phi_3(z) = z(2z-1) = 2z^2 - z
$

=== Matrices Elementales
La construcción del sistema global se basa en el cálculo previo de las matrices locales sobre el intervalo de referencia.

*Matriz de Masa ($M^e$):*
Calculando las integrales $integral_0^1 phi_i phi_j$, obtenemos:

$
  M^e = h/30 mat(
    4, 2, -1;
    2, 16, 2;
    -1, 2, 4
  )
$

*Matriz de Rigidez ($K^e$):*
Calculando las integrales de las derivadas $integral_0^1 phi_i ' phi_j '$, obtenemos:

$
  K^e = 1/(3 h) mat(
    7, -8, 1;
    -8, 16, -8;
    1, -8, 7
  )
$

=== Ensamblaje Global y Lógica de Índices

Las matrices globales se construyen mediante el proceso de ensamblaje, sumando las contribuciones locales calculadas en cada elemento. Aunque el sistema a resolver es global, el cómputo de las integrales se realiza localmente y sus valores se acumulan en las posiciones correspondientes de la matriz del sistema.

La correspondencia entre la numeración local del elemento $k$ (donde $k=1,...,N$) y la global es:
- Nodo local 1 (Izquierda) $->$ Índice global $2k - 1$
- Nodo local 2 (Centro) $->$ Índice global $2k$
- Nodo local 3 (Derecha) $->$ Índice global $2k + 1$

Esto genera una estructura donde:
1. Los nodos pares (centros) solo reciben contribución de su propio elemento (coeficiente central de las matrices locales).
2. Los nodos impares (fronteras) suman las contribuciones del elemento anterior (nodo local 3) y del actual (nodo local 1).

#pagebreak()

=== Vector de Carga ($F$)
El vector de carga se calcula elemento a elemento, transformando la integral al intervalo de referencia $[0,1]$ mediante el cambio de variable $x = x_(2k-1) + z h$.

Para cada elemento $k$, se calculan tres integrales locales numéricamente (usando Simpson, Punto Medio, etc.):
$
  I_1 = integral_0^1 f(x_(2k-1) + z h) dot phi_1(z) dot h d z \
  I_2 = integral_0^1 f(x_(2k-1) + z h) dot phi_2(z) dot h d z \
  I_3 = integral_0^1 f(x_(2k-1) + z h) dot phi_3(z) dot h d z
$

Estas contribuciones se acumulan en el vector global $F$ en las posiciones correspondientes:
$
  F(2k-1) & arrow.l F(2k-1) + I_1 \
    F(2k) & arrow.l F(2k) + I_2 \
  F(2k+1) & arrow.l F(2k+1) + I_3
$

Finalmente, se imponen las condiciones de contorno de Dirichlet nulas eliminando la primera y última fila/columna del sistema ($u_1 = u_(2N+1) = 0$).

#pagebreak()

== Resultados y Análisis de la Convergencia

Se ha implementado el método en MATLAB para mallas uniformes. Un aspecto crítico observado durante la experimentación es la influencia de la integración numérica en el cálculo del vector de carga $F$.

Para evaluar la convergencia, medimos el error en la norma $L^2(Omega)$, definido como:
$ ||u - u_h||_(L^2) = ( integral_0^1 |u(x) - u_h (x)|^2 d x )^(1/2) $
Dado que $u_h$ es un polinomio a trozos, esta integral se calcula como la suma de las integrales sobre cada elemento.

=== Comparación de las Cuadraturas

Se probaron tres métodos de integración numérica para calcular el vector de carga $F$: Rectángulo, Punto Medio y Simpson. A continuación se muestran las pendientes de convergencia ($p$) obtenidas experimentalmente ($E approx C h^p$) para cada caso:

#figure(
  table(
    columns: (auto, auto, auto),
    inset: 10pt,
    align: center,
    [*Método de Integración*], [*Pendiente Lineal ($P_1$)*], [*Pendiente Cuadrática ($P_2$)*],
    [Rectángulo (Izq)], [2.00], [2.00],
    [Punto Medio], [2.00], [2.01],
    [Simpson], [2.00], [3.00],
  ),
  caption: [Orden de convergencia observado según la regla de cuadratura.],
)

A continuación se detallan los errores absolutos específicos utilizando la cuadratura adecuada (Simpson), confirmando los órdenes óptimos:

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
  caption: [Evolución del error $L^2$ usando Simpson.],
)

==== Análisis de la degradación del orden

El error numérico total tiene dos fuentes principales: el error de aproximación del método ($E_("aprox")$) y el error cometido al calcular las integrales ($E_("int")$).
$ "Error Total" approx E_("aprox") + E_("int") $

1. *Elementos Lineales ($P_1$):* El método tiene un orden teórico $O(h^2)$. Como todas las cuadraturas usadas tienen al menos precisión $O(h^2)$, el error de integración no empeora el resultado.

2. *Elementos Cuadráticos ($P_2$):* El método debería converger a $O(h^3)$.
  - Al usar *Rectángulo o Punto Medio*, el error de integración es $O(h^2)$. Para $h$ pequeño, este error es mucho mayor que el del método ($h^2 >> h^3$), por lo que "domina" y limita la convergencia global a 2.
  - Al usar *Simpson*, el error de integración es $O(h^4)$. Al ser mucho más pequeño que el error del método ($h^4 << h^3$), se vuelve despreciable. Esto permite que se manifieste el orden real del elemento finito, recuperando la convergencia teórica $O(h^3)$.
