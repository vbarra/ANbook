# Exercices

## Exercice 1
Trouver la meilleure approximation, $\bar x$, au sens des moindres carrés, du
système
$3x=10,\qquad 4x= 5.$

Déterminer le carré de l'erreur $e^2$ et montrer que le vecteur erreur
$(10-3\bar x\;  5-4\bar x)^\top$ est orthogonal à $(3\; 4)^\top$



## Exercice 2
Trouver la droite $f(t)=at+b$ qui approche le mieux (au sens des moindres
carrés) les mesures : $f(0)=0,\quad f(1)=1,\quad f(3)=2,\quad f(4)=5$.



## Exercice 3
Soit $P=\frac{1}{\parallel v\parallel^2}vv^\top$ la matrice de 
projection sur la droite de direction $v$.

Interpréter géométriquement la matrice $Q=I-2P$. Montrer que $Q^2=I$ et
$Q^T =Q$ (\textit{i.e.} $Q$ est orthogonale). Soit $y=(y_1, y_2)^\top \in\mathbb R^2$.
Déterminer $v$ pour que la seconde coordonnée de $z=Qy$ soit nulle.


## Exercice 4
Projeter le vecteur $b=(0\; 3\; 0)^\top$ sur les droites de directions 
respectives 
$a^1=\left(\begin{array}{r}2/3\\ 2/3\\ -1/3\end{array}\right),\quad
a^2=\left(\begin{array}{r}-1/3\\ 2/3\\ 2/3\end{array}\right)$

Trouver la projection de $b$ sur le plan engendré par $a^1$ et $a^2$.


## Exercice 5
 On se donne
$A_1=\left[\begin{array}{cc}1 & -1\\ 1 & 0\\ 1 & 1\\ 1 & 2\end{array}
  \right]\quad\mbox{et }\ 
b=\left[\begin{array}{c}1 \\ 1\\ 0\\ 1 \end{array}
  \right]$

1. Mettre $A_1$ sous la forme $Q_1R_1$ où $Q_1$ est une matrice $4\times 2$ de colonnes
orthonormées $q_1$ et $q_2$, et $R_1$ une matrice $2\times 2$ triangulaire supérieure.

2.  Calculer la solution aux moindres carrés du système incompatible $A_1x_1=b$
et l'erreur $e_1$ correspondante.

3. On rajoute la colonne $a_3=[1\ 0\ 1\ 4]^T $ à la matrice $A_1$ pour avoir
$A_2=[A_1\ a_3]$. Mettre $A_2$ sous la forme $Q_2R_2$ et calculer la solution aux moindres
carrés du système $A_2x_2=b$ . 
Calculer l'erreur $e_2$ correspondante et montrer que $e_2=e_1^2-(q_3^\top b)^2$,
où $q_3$ est la troisième colonne de $Q_2$. Justifier cette diminution dans le cas général.

## Exercice 6
 On se propose de déterminer les coefficients du polynôme $p(x)=\sum_{k=0}^na_kx^k$
qui approche <<au mieux>> une fonction $y=f(x)$, connue à travers $m+1$ 
couples $(x_i,y_i)$, $i=0,\ldots,m$; avec $m\ge n$. On supposera que tous les
points de mesure $x_i$, $i=0,\ldots,m$, sont {\it distincts}. 
Pour la suite on pose
${\bf x}=(x_0,\ldots,x_m)$ et ${\bf y}=(y_0,\ldots,y_m)$.

Pour ${\bf a}=(a_0,a_1,\ldots,a_n)\in \mathbb R^{n+1}$, on définit
$J({\bf a})=\sum_{i=0}^m\Bigl[y_i-p(x_i)\Bigr]^2.$

1. Caractériser la solution ${\bf a}^\ast$ qui minimise $J({\bf a})$
sur $\mathbb R^{n+1}$. Quelle est l'interprétation du polynôme $p(x)$
quand $m=n$. Quelle est alors la valeur de $J({\bf a}^\ast)$.

2. Montrer qu'il existe une matrice $U$, rectangulaire à $(m+1)$ lignes
et $(n+1)$ colonnes, telle que, si ${\bf a}^\ast$ est un point stationnaire de
$J$, alors $U^{T}U{\bf a}^\ast=U^{T}{\bf y}$.
3. On suppose que $n=1$, i.e. $p(x)=a_0+a_1x$. On pose
$\begin{eqnarray*}
&&\bar{{\bf x}}=\frac{1}{m+1}\sum_{i=0}^mx_i,\qquad
           \bar{{\bf y}}=\frac{1}{m+1}\sum_{i=0}^my_i,\\[3pt]
&&\sigma({\bf x}^2)=\frac{1}{m+1}\sum_{i=0}^m(x_i-\bar{{\bf x}})^2,\qquad
  \sigma({\bf x},{\bf y})=\frac{1}{m+1}\sum_{i=0}^m
         \left(x_i-\bar{{\bf x}}\right)\left(y_i-\bar{{\bf y}}\right).
\end{eqnarray*}$

Montrer que
$a_1=\frac{\sigma({\bf x},{\bf y})}{\sigma({\bf x}^2)},\qquad a_0=\bar{{\bf y}}-a_1\bar{{\bf x}}$
et en déduire que le <<point moyen>> $(\bar{{\bf x}},\bar{{\bf y}})$ 
appartient à la droite $y=a_0+a_1x$.

## Exercice 7
 On se donne 
$$
A = \left[\begin{array}{ccc}1 & 2 & -1 \\ 2 & 1 &  4 \\ 1 & 1 &  1\end{array}\right]
$$ 
et 
$$
b = \left[\begin{array}{c}1 \\ 1 \\1 \end{array}\right].
$$


1. Montrer que le système $Ax = b$ est incompatible.
2. Donner une CNS pour que $x \in \mathbb R^3$
réalise le minimum de $E(x) = \| Ax - b \|_2^2$.
3. En déduire que la solution des moindres carrés n'est pas unique.  
4. Construire une base orthonormée de $ImA$ par le procédé
d'orthonormalisation de Gram-Schmidt. En déduire l'existence
d'une matrice $Q$, de format $3\times 2$, formée de vecteurs colonnes
orthonormés, et d'une matrice $R$ triangulaire supérieure de format
$2\times 3$, telles que $A = QR$.
5. Montrer que $RR^T $ est une matrice inversible. En déduire
une expression explicite pour l'ensemble des solutions des moindres
carrés de la forme $x = x_0 + \alpha u$, où $\alpha \in \mathbb R$, et
tel que $x_0$ et $u$ soient des vecteurs orthogonaux.
% \item Exprimer la solution des moindres carrés de norme minimale en fonction
% de $Q$, $R$ et $b$. On ne demande pas de la calculer mais de raisonner avec les projecteurs.
6. Exprimer $\min\| Ax - b \|_2^2$ en fonction de $b$ et d'un 
vecteur orthogonal à $\mathrm{Im}Q$.  
7. Application numérique : calculer la solution des moindres carrés de norme minimale, ainsi que l'erreur résiduelle $\min\| Ax - b \|_2^2$. 

 