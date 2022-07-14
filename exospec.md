---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
# Exercices

## Exercice 1

Soit les matrices 
$\begin{pmatrix}
1 & 1 \\
0 & 2 
\end{pmatrix}
$,
$\begin{pmatrix} 
1 & 0 \\
1 & 2 
\end{pmatrix}
$,
$\begin{pmatrix} 
1 & 1 \\
0 & 1 
\end{pmatrix}
$,
$
\begin{pmatrix}
1 & 1 \\
1 & 0 
\end{pmatrix}
$,
$
\begin{pmatrix}
1 & 2 \\
2 & 4 
\end{pmatrix}
$,
$
\begin{pmatrix}
0 & 1 & 0 \\
1 & 0 & 1 \\
0 & 1 & 0
\end{pmatrix}
$,
$
\begin{pmatrix} 
11 & -5 & 5 \\
-5 & 3 & -3 \\
5 & -3 & 3
\end{pmatrix}
$

Pour chacune de ces matrices, répondre aux questions suivantes :
- Polynôme caractéristique ?
- Spectre et rayon spectral ?
- Multiplicité algébrique des valeurs propres ?
- Sous-espaces propres ?
- Multiplicité géométrique des valeurs propres ?
- La matrice est-elle diagonalisable ?  Si oui, la diagonaliser.


## Exercice 2
Que peut-on dire de la diagonalisabilité, des valeurs propres et des vecteurs propres
- d'une homothétie ?
- d'un projecteur ?
- d'une rotation?


## Exercice 3

1. Soit $A=\begin{pmatrix}  
1 & 1 \\
0 & 2 
\end{pmatrix}$.
Exprimer  la matrice $A^k$.

2.  Soit la suite de Fibonacci $1$, $1$, $2$, $3$, $5$, $8$, $13$, ... définie par :
$ \begin{align*}
u_0&=u_1=1,\\
u_{k}&=u_{k-1}+u_{k-2}, \quad k\geq 2.
\end{align*}$

En posant $X_k=\begin{pmatrix}  u_{k+1} \\ u_k \end{pmatrix}^\top $, déterminer la matrice $A$ telle que $X_k=A X_{k-1} $ puis en déduire la limite de $F_k=\frac{u_{k+1}}{u_k}$ quand $k$ tend vers l'infini.


3. Soientt les suites récurrentes linéaires simultanées suivantes :
$
\left\lbrace\begin{array}{l}
u_0=v_0=0,w_0=1\\
u_{n+1}=11u_n-5v_n+5w_n\\
v_{n+1}=-5u_n+3v_n-3w_n\\
w_{n+1}=5u_n-3v_n+3w_n
\end{array}\right.
$

Exprimer $u_n,v_n,w_n$ en fonction de $n$, pour $n\geq 0$. 

4.  Résoudre le système différentiel suivant (variable $t$, inconnues $x(t), y(t), z(t)$ à valeurs réelles)
$
\left\lbrace\begin{array}{l}
x'(t)=y(t)\\
y'(t)=x(t)+z(t)\\
z'(t)=y(t)
\end{array}\right.
$
}


## Exercice 4

1. Soit $a$ et $b$ deux réels et $D$ une matrice carrée. Montrer que si
$\lambda$ est valeur propre de $D$, alors $a\lambda+b$ est une valeur
propre de $aD+b\mbb I$, o\`u $\mbb I$ est la matrice identité.

2. Soit $D$ la matrice carrée tridiagonale définie par
$\begin{eqnarray*}
&&d_{ii}=0,\quad i=1,\ldots,n\\
&&d_{i,i+1}=d_{i+1,i}=1,\quad i=1,\ldots,n-1\\
&&d_{ij}=0,\quad\mbox{ailleurs}.
\end{eqnarray*}$

Pour $j=1,\ldots,n$, soit $\mbf x^{(j)}$ le vecteur de $\mbb R^n$ dont
la $i$--ème composante est
$
x_i^{(j)}=\sin 2i\frac{j\pi}{n+1}.
$

Montrer que $\mbf x^{(j)}$ est un vecteur propre de $D$ associée à la
valeur propre 
$
\lambda_j=2\cos 2\frac{j\pi}{n+1}.
$

3. En déduire les vecteurs propres et les valeurs propres de la
matrice $4\times 4$ suivante
$
A=\left[\begin{array}{rrrr}
  2 & -1 &  0 &  0\\
       -1 &  2 & -1 &  0\\
        0 & -1 &  2 & -1\\
        0 &  0 & -1 &  2
        \end{array}
  \right]
$

 

 

## Exercice 5
Soit $a$ et $b$ deux vecteurs non colinéaires de $\mbb R^n$, $\norme{a}=\norme{b}=1$. Déterminer
les vecteurs propres et les valeurs propres de la matrice $n\times n$
$
A=aa^\top +bb^\top .
$

 

## Exercice 6
Soit $A'=A+E$, la somme de deux matrices symétriques. Les valeurs propres
de ces trois matrices sont notées
$
\lambda'_1\le \lambda'_2\le\cdots\le\lambda'_n,\quad
\lambda_1\le \lambda_2\le\cdots\le\lambda_n,\quad
\mu_1\le\mu_2\le\cdots\le\mu_n.
$

1. Montrer les propriétés suivantes pour tout $i=1,\ldots,n$
- $\lambda_i+\mu_1\le\lambda'_i\le\lambda_i+\mu_n$
- $|\lambda'_i-\lambda_i|\le \parallel E\parallel$, quelle que soit 
la norme matricielle $\parallel\cdot\parallel$.


2. Soit $E^{(k)}=A^{(k)}-\mathrm{diag}(a_{ii}^{(k)})$, où les
$A^{(k)}$ sont les matrices engendrées par la méthode de Jacobi.
En utilisant la norme $\parallel E\parallel_F=(\mbox{trace}(E^\top E))^{1/2}$,
montrer que $E^{(k)}$ tend vers la matrice nulle, lorsque 
$k\longmapsto \infty$.

3. En déduire le théorème de convergence de la méthode de Jacobi, i.e.
- $a_{ij}^{(k)}\longrightarrow 0$, pour $i\ne j$
- chaque $a_{ii}^{(k)}\longrightarrow \lambda_i$, o\`u $\lambda_i$
est une valeur propre de $A$. 


## Exercice 7
Supposons calculée la plus grande valeur propre $\lambda_n$ d'une 
matrice symétrique $A$ et le vecteur propre associé $v_n$.

1. Quelle est la matrice $P_n$ de projection orthogonale sur 
$v_n^\perp$.

2. Montrer que $P_nA$ a les m\^emes vecteurs propres que $A$, à
l'exception de $v_n$. Quelles sont les valeurs propres de $P_nA$.
En déduire une première méthode de calcul des plus grandes valeurs 
propres de $A$.
\end{exo}


## Exercice 8
 Chaque année, une proportion $a$ de clermontois quitte Clermont-Ferrand pour s'installer 
dans le Puy de Dôme (hors Clermont-Ferrand), une proportion $b$ de non clermontois
du Puy de Dôme s'installe à Clermont-Ferrand ($0\le a,b\le 1$). On pose $x_0$ le nombre
de Clermontois, et $y_0$ le nombre de non clermontois (du Puy de Dôme), en cette année.
 Au bout de $k$ années, ces nombres vaudront $x_k$ et $y_k$. On suppose qu'il n'y a pas
de flux migratoire du Puy de Dôme vers l'extérieur, et vice versa.

1. Calculer $x_k$ et $y_k$ en fonction de $a$, $b$, $x_0$ et $y_0$.
2. A quelle condition sur $a$ et $b$ existe-t-il une limite pour $x_k$ et
$y_k$ lorsque $k\rightarrow+\infty$? Que vaut cette limite?
3. Application numérique: calculer cette limite en fonction de $x_0$ et $y_0$ lorsque
$a=0.1$ et $b=0.2$.


