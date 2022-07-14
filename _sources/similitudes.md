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
# Similitudes

L'objectif est encore une fois de transformer une matrice par des transformations simples en une matrice dont on connaît les valeurs propres, c'est-à-dire, une matrice triangulaire ou diagonale. Les transformations qui maintiennent le spectre d'une matrice sont des similitudes.


```{prf:definition} Similitude
Deux matrices carrées $A$ et $B$ sont dites *semblables* 
s'il existe une matrice $S$ non singulière telle que 

$B=S^{-1}AS$
```
```{index} Matrice;semblable
```

La transformation de $A$ vers $B$ est une *similitude*. En l'écrivant sous la forme $AS=SB$, on retrouve une généralisation de la définition des valeurs propres et des vecteurs propres. On a d'ailleurs le résultat fondamental :

```{prf:property}
Deux matrices semblables ont les mêmes valeurs propres
```
En effet, soit $x$ un vecteur propre de $A$ associé à la valeur propre $\lambda$. On a donc $Ax=\lambda x$, qui s'écrit $SBS^{-1}x=\lambda x$, ce qui veut dire que $\lambda$ est valeur propre de $B$ associé au vecteur propre $S^{-1}x$.



```{code-cell} ipython3
from sympy import Matrix
A = Matrix([[2, 1], [1, 2]])
S = Matrix([[1, 4], [0, 1]])
```
Pour $\bf S$ inversible, on construit donc $\bf{B=S^{-1}AS}$ semblable. 

```{code-cell} ipython3
B = S.inv() * A * S
```

On vérifie que $\bf A$ et $\bf B$ ont le même spectre
```{code-cell} ipython3
A.eigenvals(), B.eigenvals()
```

Les valeurs propres sont donc 1 et 3, de multiplicité 1 chacune.

En revanche, les vecteurs propres sont différents

```{code-cell} ipython3
Matrix(A.eigenvects())
Matrix(B.eigenvects())
```



L'intérêt de ces transformations est double : 
- les valeurs propres sont inchangées
- en supposant les vecteurs propres linéairement indépendants, la similitude associée à la matrice $X$ dont les colonnes sont les vecteurs propres transforme $A$ en une matrice diagonale dont les éléments diagonaux sont les valeurs propres de $A$ : $X^{-1}AX = \Lambda$.


Montrons maintenant sur un exemple que deux matrices semblables représentent la même transformation linéaire sur deux bases différentes : soit $P$ la matrice de projection dans $\mathbb{R}^2$ sur la droite $L$ d'angle $\theta$ : 
$P =
\left[
\begin{array}{cc}
\cos^2(\theta) & \cos(\theta)\sin(\theta)\\
\cos(\theta)\sin(\theta)& \sin^2(\theta)
\end{array}
\right]
=uu^\top \quad\mbox{avec}\quad u=\begin{pmatrix}\cos(\theta)\\ \sin(\theta)\end{pmatrix}
$

Si on fait maintenant tourner la base canonique orthonormée d'un angle $\theta$, la projection devient maintenant une projection sur l'axe horizontal et s'écrit

$Q =
\left[
\begin{array}{cc}
1&0\\
0&0\\
\end{array}
\right]
$

Pour passer de $P$ à $Q$, on utilise la matrice de rotation d'angle $\theta$ 

$R =
\left[
\begin{array}{cc}
\cos(\theta) & -\sin(\theta)\\
\sin(\theta)& \cos(\theta)
\end{array}
\right]
$

Le changement de base se traduit par : si $x$ est le vecteur de coordonnées dans la base de départ, et $u$ le vecteur de coordonnées dans la base d'arrivée, on a 
$=Ru$.

Soit $y$ la projection de $x$ sur $L$ dans la base de départ, et $v$ ses coordonnées dans la base d'arrivée. Alors 

$\begin{align*}
y &= Px\\
Rv&= PRu\\ \Rightarrow v&=R^{-1}PRu=Qu
\end{align*}$

et donc $Q=R^{-1}PR$ et $P$ et $Q$ ont les mêmes valeurs propres. Comme $Q$ est diagonale, on retrouve ici le résultat connu : toute matrice de projection sur une droite dans $\mathbb{R}^2$ a pour valeurs propres 1 et 0, et pour vecteurs propres associés les colonnes de la matrice de rotation $R$, c'est-à-dire les directions de $L$ et $L^{\bot}$ respectivement.