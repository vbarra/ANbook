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
# Méthode de Jacobi

L'intérêt principal des matrices réelles symétriques est qu'il existe une base de vecteurs propres orthonormés. On peut donc la diagonaliser par une transformation orthogonale : soient $A$ une telle matrice et $Q$ la matrice orthogonale dont les colonnes sont les vecteurs propres de $A$, alors $Q^\top AQ=diag\{\lambda_1\cdots \lambda_n\}$.

La méthode de Jacobi est une méthode d'élimination symétrique itérative utilisant des similitudes
- la matrice transformée tend vers une matrice diagonale
- le produit des transformations orthogonales tend vers la matrice des vecteurs propres.


## Principe d'élimination symétrique

Ce principe, basé sur des rotations successives, sera illustré tout d'abord sur une matrice (2$\times$2) : soit $A$ la sous-matrice (2$\times$2) symétrique 

$A=\begin{pmatrix}a_{pp}\quad a_{pq}\\
         a_{pq}\quad a_{qq}
   \end{pmatrix}\quad a_{pq}\neq 0$

Soit $R$ la matrice de rotation d'angle $-\theta$ : 

$R = \begin{pmatrix}\cos\theta \quad  \sin\theta\\ 
                  -\sin\theta\quad  \cos\theta
    \end{pmatrix}$

Pour éliminer l'élément $a_{pq}$ non diagonal, on détermine $\theta$ tel que 

$R^\top AR=\begin{pmatrix}*\quad0\\0\quad *\end{pmatrix}$

Le calcul donne 
$\mathrm{cotg}(2\theta)=\frac{a_{qq}-a_{pp}}{2a_{pq}}$

Considérons maintenant une matrice symétrique $n\times n$ $A$ telle que l'élément $a_{pq}$ soit non nul. La transformation orthogonale $\Omega$ telle que $A'=\Omega^\top A\Omega$, avec $a'_{pq}=a'_{qp}=0$ est :

$\Omega =
\left (
\begin{array}{ccccc}
1&&&&\\
&\ddots& \cos\theta & \sin\theta&\\
&&-\sin\theta & \cos\theta&\\
&&&\ddots&\\
&&&&1\\
\end{array}
\right )
\left .
\begin{array}{c}
p\\
q\\
\\
\end{array}
\right .$

On vérifie que seules les lignes et colonnes $p$ et $q$ de $A$ sont modifiées. 
En posant $c=\cos\theta$, $s=\sin\theta$, $t=\tan\theta$, la mise à jour s'écrit :

$\begin{align*}
a'_{pj}&=ca_{pj}-sa_{qj},j\neq p,q\\
a'_{qj}&=ca_{qj}+sa_{pj},j\neq p,q\\
a'_{pp}&=a_{pp}-ta_{pq}\\
a'_{qq}&=a_{qq}+ta_{pq}
\end{align*}$

Un choix classique pour $p$ et $q$ est celui qui permet d'éliminer l'élément non diagonal de plus grand module : $|a_{pq}|=\displaystyle\max_{i\neq j}|a_{ij}|$


## Convergence

Soit $\{A_k\}$ la suite de matrices engendrée par l'algorithme en éliminant à chaque itération un élément non diagonal non nul ($A_0=A$) et son symétrique. Comme la transformation est orthogonale, la norme de Frobenius de la matrice est conservée : 

$\displaystyle\sum{i,j}a_{ij}'^2=\displaystyle\sum_{i,j}a_{ij}^2$

Mais ce même résultat est vrai pour la matrice (2$\times$2) transformée ci-dessus. Donc :

$a'^2_{pp}+a'^2_{qq}=a^2_{pp}+a^2_{qq}+2a_{pq}^2$

Comme seules les lignes $p$ et $q$ sont modifiées, ces résultats impliquent que la somme des carrés des éléments diagonaux augmente strictement de la valeur $2a_{pq}^2$, et que la somme des carrés des éléments non diagonaux diminue strictement de la même quantité.

On démontre alors le théorème suivant :
```{prf:theorem} Convergence de la méthode de Jacobi
La suite de matrices $\{A_k\}$ engendrée par la méthode de Jacobi classique converge vers une matrice diagonale contenant toutes les valeurs propres de $A$ sur la diagonale
```


On peut également démontrer que, si les valeurs propres sont toutes distinctes, la suite des matrices $Q_k=\Omega_1\cdots \Omega_k$ converge vers la matrice orthogonale contenant les vecteurs propres.

```{warning}
Un élément annulé peut redevenir non nul aux itérations suivantes.
```

```{prf:remark} 
:class: dropdown
- Le tri du plus grand élément parmi $n(n-1)/2$ coûtant relativement cher, on lui préfère d'autres stratégies plus économiques (balayage cyclique ou choix avec seuil)
- La méthode de Jacobi présente d'excellentes performances pour des matrices pleines de faible dimension (typiquement inférieure à 100). Dans le cas général, on lui préférera la méthode QR.
``` 



