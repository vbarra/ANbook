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


# Gauss et LU
Parmi toutes les méthodes de résolution de systèmes linéaires, nous étudions dans la suite deux algorithmes largement utilisés : l'élimination de Gauss et la méthode LU.

# Préambule : transformations élémentaires
## Définition
````{prf:definition} Transformation élémentaire
```{index} Transformation ; élémentaire
```
Soit ${\bf A}\in\mathcal{M}_{mn}(\mathbb R)$. On appelle transformation élémentaire
des lignes (respectivement des colonnes) de ${\bf A}$ une combinaison des trois transformations
suivantes :
- multiplication d'une ligne (resp. colonne) par un scalaire non nul;
- remplacement d'une ligne(resp. colonne)  par la somme de cette ligne et d'une autre ligne
  (resp. colonne);
- permutation de deux lignes (resp. colonnes).
````

Une transformation élémentaire sur les lignes (resp. colonnes) de ${\bf A}$
équivaut à multiplier ${\bf A}$ à gauche (resp. à droite) par une matrice
élémentaire obtenue en appliquant la même transformation à la matrice
identité (de taille $m$ pour les transformations sur les lignes, et de taille $n$ pour les transformations sur les colonnes). C'est une transformation non singulière, elle ne change donc pas le rang de ${\bf A}$.


Si on note pour $\alpha\in\mathbb R^*$:

- ${\bf M_i(\alpha)}$ la matrice élémentaire qui réalise ${\bf A_{i,\bullet}}\leftarrow \alpha {\bf A_{i,\bullet}}$ (ou ${\bf A_{\bullet,i}}\leftarrow \alpha {\bf A_{\bullet,i}}$)
- ${\bf C_{ij}(\alpha)}$ la matrice élémentaire qui réalise ${\bf A_{i,\bullet}}\leftarrow {\bf A_{i,\bullet}} + \alpha {\bf A_{j,\bullet}}$ (ou ${\bf A_{\bullet,i}}\leftarrow {\bf A_{\bullet,i}} + \alpha {\bf A_{\bullet,j}})$
- ${\bf P_{ij}}$ la matrice élémentaire qui réalise ${\bf A_{i,\bullet}}\leftrightarrow {\bf A_{j,\bullet}}$  (ou ${\bf A_{\bullet,i}}\leftrightarrow {\bf A_{\bullet,j}}$)

Alors on a le résultat : 


```` {prf:theorem} 
Les matrices ${\bf M_i(\alpha)}$, ${\bf C_{ij}(\alpha)}$ et ${\bf P_{ij}}$ sont inversibles et :

$({\bf M_i(\alpha))^{-1}} = {\bf M_i(\alpha^{-1})}\quad ({\bf C_{ij}(\alpha))^{-1}} =  {\bf C_{ij}(-\alpha)} \quad {\bf P_{ij}^{-1}}=  {\bf P_{ij}}$
````

````{prf:example}
Soit $
{\bf A}=
\begin{pmatrix}1 & 0 & 0 & 2\\ 1 & 2 & 0 & 1\\ 2 & 1 & 3 & 4\end{pmatrix}
$

Si on additionne les deux premières lignes et que l'on inscrive le résultat
sur la première ligne sans changer les deux autres, on obtient
$
{\bf A'}=
\begin{pmatrix}2 & 2 & 0 & 3\\ 1 & 2 & 0 & 1\\ 2 & 1 & 3 & 4\end{pmatrix}
$

On peut vérifier que ${\bf A'}={\bf C_{12}(1)}{\bf A}$, avec
$
{\bf C_{12}(1)}=\begin{pmatrix}1 & 1 & 0 \\ 0 & 1 & 0\\ 0 & 0 & 1\end{pmatrix}
$

De même la permutation des première et dernière colonnes de ${\bf A}$ s'effectue en multipliant à droite ${\bf A}$ par 
${\bf P_{14}}=
\begin{pmatrix}0&0&0 & 1  \\ 0 & 1 & 0&0\\ 0 & 0 & 1&0\\1&0&0&0\end{pmatrix}
$
````

## Pivotage d'une colonne
```{index} Pivotage
```

C'est une séquence de transformations élémentaires qui transforme une colonne
d'une matrice ${\bf A}$ en vecteur unitaire. L'élément qui doit se transformer en 1 est le `pivot`.
```{index} Pivot
```
L'algorithme {prf:ref}`pivotage` décrit la transformation d'une colonne ${\bf A_{\bullet s}}$ en le $s^e$ vecteur de la base canonique.

```{prf:algorithm} Pivotage d'une colonne ${\bf A_{\bullet s}}$ suivant le pivot $a_{rs}\ne 0$
:label: pivotage

**Entrée : ** ${\bf A}$, $s$ l'indice de la colonne à pivoter, $r$ l'indice de ligne du pivot

**Sortie : ** La colonne transformée

1. ${\bf A_{r \bullet}}\leftarrow \frac{1}{a_{rs}}\bf A_{r \bullet}$
2. Pour $i \ne r$ : ${\bf A_{i \bullet}}\leftarrow {\bf A_{i \bullet}} -a_{is}{\bf A_{r \bullet}}$
```



Pour $m=n$, il est clair que, si tous les pivots successifs sont non nuls, une séquence de
$n$ pivotages effectués sur les $n$ colonnes de la matrice avec le pivot
sur la diagonale transformera la matrice en l'identité. Si ${\bf E_i}$ est la 
matrice élémentaire associée au pivotage de la $i^e$ colonne, on a
$
\mathbb I_n={\bf E_n}{\bf E_{n-1}}\cdots {\bf E_1}{\bf A}
$
d'où une première méthode pour calculer l'inverse d'une matrice : 
$
{\bf A^{-1}}={\bf E_nE_{n-1}}\cdots {\bf E_1}.
$

Donc pour calculer l'inverse d'une matrice, il suffit d'effectuer les pivotages
en parallèle sur la matrice identité. La matrice à inverser se transforme
progressivement en la matrice identité alors que l'identité devient 
l'inverse. 

S'il est vrai que, quand aucun pivot nul n'est rencontré, cet
algorithme a une complexité de $O(n^3)$, il peut atteindre $O(2n^3)$ dans le
cas général et nécessite le stockage de deux matrices $n\times n$. On lui
préférera les méthodes étudiées dans la suite, basées elles aussi sur
des pivotages successifs de la matrice, mais plus robustes et moins 
coûteuses.

````{prf:example}
Pivotons la première colonne de 
$
{\bf A}=
\begin{pmatrix}1 & 1 & 1\\ 1 & 0 & -1\\ -2 & 1 & 3\end{pmatrix}
$
avec $a_{11}=1$ comme pivot. On doit donc multiplier $A$ à gauche par
la matrice élémentaire 
$
{\bf E}=
\begin{pmatrix}1 & 0 & 0\\ -1 & 1 & 0\\ 2 & 0 & 1\end{pmatrix}$
et le résultat est 
${\bf A'}={\bf E}{\bf A}=
\begin{pmatrix}1 & 1 & 1\\ 0 & -1 & -2\\ 0 & 3 & 5\end{pmatrix}
$
````

# Elimination de Gauss
## Exemple introductif


 Soit le système de 3 équations à 3 inconnues
\begin{eqnarray}
2x_1-3x_2\phantom{-3x_3}&=& 3\label{SL1}\\ 
4x_1-5x_2+x_3& =& 7\label{SL2}\\ 
2x_1-x_2-3x_3& =& 5\label{SL3}
\end{eqnarray}
La méthode d'élimination de Gauss (ou méthode du pivot) consiste à utiliser la
première équation pour calculer $x_1$ en fonction des autres variables puis de
remplacer cette variable dans les équations suivantes. Cette élimination se
poursuit avec $x_2$ dans les nouvelles équations (sauf la première) jusqu'à
l'obtention d'une équation à une seule inconnue. On remonte alors en remplaçant
les variables calculées dans les équations ayant servi à l'élimination :

Tirons $x_1$ de l'équation (\ref{SL1}):
\begin{equation}x_1=\frac{3}{2}(1+x_2) \label{VarX}\end{equation}
Remplaçons $x_1$ dans les deux dernières équations (\ref{SL2})-(\ref{SL3})
par son expression (\ref{VarX})
\begin{eqnarray}
x_2+x_3 &=& 1 \label{SLR1}\\ 
2x_2-3x_3 &=& 2\label{SLR2}.
\end{eqnarray}
Tirons $x_2$ de l'équation (\ref{SLR1}) de ce nouveau système :
\begin{equation}x_2=1-x_3. \label{VarY}\end{equation}
 Remplaçons $x_2$ dans l'équation (\ref{SLR2}) :
 \begin{equation}-5x_3=0.\label{VarZ}\end{equation}
 La phase d'élimination est terminée. On effectue alors la substitution en sens
inverse des variables gr\^ace aux équations (\ref{VarZ}), (\ref{VarY}) et
(\ref{VarX}):
\begin{eqnarray*}
(\ref{VarZ}) &\Rightarrow& x_3=0\\ (\ref{VarY}) &\Rightarrow& x_2=1\\ (\ref{VarX})
&\Rightarrow& x_1=3
\end{eqnarray*}

On remarque qu'à une étape donnée, l'élimination d'une variable peut se faire
dans n'importe quelle équation, à condition, bien sûr, que cette équation
contienne la variable en question.


## Systèmes triangulaires

On constatera au paragraphe suivant que la phase d'élimination consiste à
transformer le système original en un système triangulaire. Un système
triangulaire est un système dont la matrice est triangulaire.

Un système triangulaire supérieur se résout par substitution inverse
\index{substitution!inverse}
(on supposera tous les éléments diagonaux $a_{ii}$ non nuls).
\begin{eqnarray*}
a_{11}x_1+a_{12}x_2+\cdots+a_{1n}x_n&=&b_1\\ a_{2n}x_2+\cdots+a_{2n}x_n &=&
b_2\\ \vdots & & \vdots\\ a_{nn}x_n&=&b_n
\end{eqnarray*}

On commence d'abord par calculer $x_n$ dans l'équation $n$. Puis à l'aide de
$x_n$, on peut calculer $x_{n-1}$ dans l'équation $n-1$ et ainsi de suite
jusqu'à $x_1$ : 
\begin{eqnarray}
x_n    &=& \frac{b_n}{a_{nn}},\label{BckSub1}\\[.35pc]
x_{k}&=&\frac{1}{a_{kk}}\left[b_k-\sum_{j=k+1}^na_{kj}x_j\right],\quad
k\in[\![n-1, 1]\!]\label{BckSub2}
\end{eqnarray}
On remarque que le calcul de $x_k$ coûte $n-k$ flops et une division. Le
coût total de l'algorithme est donc de $$1+2+\cdots+n-1=\frac{n(n-1)}{2},$$
soit (on ne garde que les termes de plus haut degré)
$$
\frac{{n^2}}{ 2}\quad\mbox{ flops}\quad\mbox{ et}\quad 
 n\quad\mbox{ divisions.}
$$

Dans le cas d'un système triangulaire inférieur, on effectue des substitutions
directes \index{substitution!directe} : 
\begin{eqnarray}
x_1    &=& \frac{b_1}{a_{11}},\label{FwdSub1}\\[.35pc]
x_{k}&=&\frac{1}{a_{kk}}\left[b_k-\sum_{j=1}^{k-1}a_{kj}x_j\right],\quad
k\in[\![2, n]\!].\label{FwdSub2}
\end{eqnarray}
Le coût de l'algorithme de substitution (\ref{FwdSub1})-(\ref{FwdSub2}) est
le même que celui de substitution inverse (\ref{BckSub1})-(\ref{BckSub2}).

## Méthode de Gauss (ou du pivot) pour les systèmes linéaires


Montrons d'abord que la technique d'élimination 
correspond à une opération de pivotage.

Eliminer $x_1$ de l'équation (\ref{SL2}) revient à faire la somme de l'équation
(\ref{SL1}) multipliée par -2 avec (\ref{SL2}). De même éliminer $x_1$ de
(\ref{SL3}) revient à faire la somme de l'équation (\ref{SL1}) multipliée par
-1 avec l'équation (\ref{SL3}). On a donc effectué un pivotage de la première
colonne (celle de $x_1$) en utilisant la première ligne comme ligne du pivot.\index{pivot}

\begin{rem}\rm Lors de l'exécution de l'algorithme : 
\begin{enumerate} 
\item A chaque itération, l'élément de la colonne dans la ligne du pivot
doit être non nul.
\item Les transformations élémentaires effectuées sur la matrice sont
effectuées en parallèle sur le second membre.
\item A la fin de l'élimination (si tout se passe bien!), on obtient un système
triangulaire avec les pivots sur la diagonale.
\end{enumerate}
\end{rem}

\exemple{
Les transformations successives sont décrites ci-dessous pour l'exemple
du  paragraphe \ref{UnExemple}.
\noindent
Itération 1 :
$$
\left(\begin{array}{rrr} 1 & & \\ -2& 1& \\ -1 & & 1 \end{array}\right)
\left(\begin{array}{rrr|r}
       (2)& -3 & 0 & 3\\ 4 & -5 & 1 & 7\\ 2 & -1 & -3 & 5\\
      \end{array} 
\right)
=\left(\begin{array}{rrr|r}
       2& -3 & 0 & 3\\ 0 & 1 & 1 & 1\\ 0 & 2 & -3 & 2\\
      \end{array} 
\right)
$$
\noindent
Itération 2 :
$$
\left(\begin{array}{rrr} 1 & & \\   & 1& \\  &  -2 & 1 \end{array}\right)
\left(\begin{array}{rrr|r}
       2 & -3 & 0 & 3\\ 0 & (1) & 1 & 1\\ 0 & 2 & -3 &  2\\
      \end{array} 
\right)
=\left(\begin{array}{rrr|r}
       2& -3 & 0 & 3\\ 0 & 1 & 1 & 1\\ 0 & 0 & -5 & 0\\
      \end{array} 
\right)
$$
}

On pose alors 
\begin{eqnarray*}
{\bf A^{(1)}}&=&{\bf A}\\[.2pc]
{\bf b^{(1)}}&=&{\bf b}.
\end{eqnarray*}
Dans le cas général, on a donc à l'étape $k$ de l'algorithme une matrice
${\bf A^{(k)}}\in\mathcal{M}_{n,n+1}(\mathbb R)$ (la dernière colonne représente
le second membre du système transformé, notée ${\bf b^{(k)}}$) dont les $k-1$ premières colonnes
sont triangulaires supérieures (cf. figure \ref{Gaussétapek}).
\bigskip
\begin{figure}[htp]
\begin{center}
\unitlength=1mm
\begin{picture}(80,60)(0,0)
\put(0,0){\line(1,0){70}}
\put(70,0){\line(0,1){60}}
\put(60,0){\line(0,1){60}}
\put(0,60){\line(1,0){70}}
\put(0,0){\line(0,1){60}}
\put(25,0){\line(0,1){30}}
\put(25,30){\line(1,0){45}}
\put(25,30){\line(-5,6){25}}
\put(25,25){\framebox(10,5){$a^{(k)}_{kk}$}}
\put(8,20){\makebox(10,10){\Large (O)}}
\put(30,35){\makebox(10,10){\Huge *}}
\put(-20,25){\makebox(10,10){$\left[{\bf A^{(k)}}|{\bf b^{(k})}\right]=$}}
\end{picture}
\caption{Etape $k$ de la méthode de Gauss} 
\label{Gaussétapek}
\end{center}
\end{figure}


On met donc à zéro chaque élément de la colonne $k$ sous le pivot\index{pivot}
$a^{(k)}_{kk}$ (supposé non nul) en remplaçant la ligne $i$ (pour $i\in[\![k+1, n]\!]$) par
%$$
%\mbox{ligne}\,i\ \leftarrow\ 
%\mbox{ligne}\,i+\frac{-a^{(k)}_{ik}}{a^{(k)}_{kk}}\times\mbox{ligne}\,k.
%$$
$$
{\bf A^{(k)}_{i \bullet}}\ \leftarrow\ 
{\bf A^{(k)}_{i \bullet}}+\frac{-a^{(k)}_{ik}}{a^{(k)}_{kk}}\times {\bf A^{(k)}_{k \bullet}}
$$
L'élément $a^{(k)}_{ij}$ pour $i\in [\![k+1,n]\!]$ et $j\in [\![k+1,n+1]\!]$ devient
$$
a^{(k+1)}_{ij}=a^{(k)}_{ij}-\frac{a^{(k)}_{ik}}{a^{(k)}_{kk}}a^{(k)}_{kj}
$$
La mise à jour de la ligne $i$ requiert donc $n-k+1$ flops et 1 division. Donc
l'étape $k$ coûte $(n-k)(n-k+1)$ flops et $n-k$ divisions. On en déduit le 
coût total de l'algorithme de Gauss:
\begin{eqnarray*}
n(n-1)+(n-1)(n-2)+\cdots+2\cdot 1&=&n^2+(n-1)^2+\cdots+2^2+1^2
                                     -(n+n-1+\cdots+2+1)\\
 &=&\frac{1}{6}n(n+1)(2n+1)-\frac{1}{2}n(n+1)\ \mbox{flops}
\end{eqnarray*}
et 
$$
\frac{1}{2}n(n+1)\ \mbox{divisions.}
$$
On dira que la complexité de la méthode de Gauss est de 
$\displaystyle{\frac{1}{3}{n^3}}$ flops.\index{flops}

L'algorithme d'élimination de Gauss est présenté dans l'algorithme~\ref{MGauss1}.
Le second membre est stocké dans la dernière colonne de ${\bf A}$. L'élimination
transforme la matrice ${\bf A}$ en une matrice triangulaire supérieure (\textit{étape 1.}).
Le système triangulaire supérieur est ensuite résolu par substitution inverse
(\textit{étape 2.}).

\begin{algorithm}
\caption{Méthode de Gauss}
%%-----------------------------------------------
\Entree{${\bf A}\in\mathcal{M}_{n,n+1}(\mathbb R)$ (la dernière colonne représente
le second membre du système)}
\Sortie{Solution ${\bf x}$ du syst\`eme ${\bf Ax}={\bf b}$}
\Deb{
\textit{étape 1. Elimination}\;
\Pour{$k=1,\ldots, n-1$}{
    \Pour{$i=k+1,\ldots, n$}{
    \Pour{$j=k+1,\ldots,n+1$}{
         $\dps{a_{ij} \leftarrow  a_{ij}-\frac{a_{ik}}{a_{kk}}a_{kj}}$\;
        }
    }
}
\textit{étape 2. Résolution du système triangulaire}\;
$\dps{x_n\leftarrow \frac{a_{n,n+1}}{a_{nn}}}$\;
\Pour{$k=n-1,\ldots,1$}{
    $\dps{x_k\leftarrow \frac{1}{a_{kk}}\left[
                                  a_{k,n+1}-\sum_{j=k+1}^na_{kj}x_j
                                               \right]}$\;
}}
\label{MGauss1}
\end{algorithm}


## Stratégies pour les pivots nuls


Si le pivot est nul (en pratique, on évitera aussi les pivots de valeur
absolue trop petite, cf. chapitre 3), on le remplace en effectuant une 
permutation avec un 
élément non nul parmi les éléments sous lui et/ou à sa droite. On distingue
généralement deux stratégies :
\begin{description}
\item[Pivotage total] \index{pivot!total}On permute lignes et colonnes pour 
choisir le plus grand élément en valeur absolue dans la sous-matrice en bas à 
droite (voir
Fig. \ref{Gaussétapek}) :
$$
\max\left\{\left|a^{(k)}_{ij}\right|\ ; \ i,j\in [\![k,n]\!]
    \right\}
$$
\item[Pivotage partiel] \index{pivot!partiel}On ne permute que les lignes sous le 
pivot en  choisissant le plus grand élément en valeur absolue
$$
\max\left\{\left|a^{(k)}_{ik}\right|\ ; \ i\in [\![k,n]\!]
    \right\}
$$
\end{description}

La stratégie du pivotage partiel est la plus utilisée car la plus économique
(la recherche du plus grand élément sur une liste de $p$ nombres coûte
$p$ comparaisons numériques, ce qui donne $n^3/3$ comparaisons pour le pivotage
total). L'algorithme d'élimination de Gauss avec recherche du pivot partiel est 
décrit dans l'algorithme~\ref{MGauss2}. Le second membre est
stocké dans la dernière colonne de ${\bf A}$.  

\begin{algorithm}
\caption{Méthode de Gauss avec pivot partiel}
\label{MGauss2}
\Entree{${\bf A}\in\mathcal{M}_{n,n+1}(\mathbb R)$ (la dernière colonne représente
le second membre du système)}
\Sortie{Solution ${\bf x}$ du syst\`eme ${\bf Ax}={\bf b}$}
\Deb{
\textit{1. Elimination}\;
\Pour{$k=1,\ldots, n-1$}{
         \textit{Recherche du pivot}\;
         $c_p\leftarrow |a_{kk}|,\: i_p\leftarrow k$\;
         \Pour{$i=k+1,\ldots,n$}{
         \Si{$ |a_{ik}|> c_p$}{
           $c_p\leftarrow |a_{ik}|,\: i_p\leftarrow i$
         }
        }
         \textit{Permutation}\;
         \Si{$i_p\ne k$}{
         Permuter les lignes $i_p$ et $k$ de la matrice $A$
         }
         \textit{Pivotage}\;
    \Pour{$i=k+1,\ldots, n$}{
    \Pour{$j=k+1,\ldots,n+1$}{
         $\dps{a_{ij} \leftarrow  a_{ij}-\frac{a_{ik}}{a_{kk}}a_{kj}}$\;
      }
    }
}
\textit{2. Résolution du système triangulaire}\;
$\dps{x_n\leftarrow \frac{a_{n,n+1}}{a_{nn}}}$\;
\Pour{$k=n-1,\ldots,1$}{
    $\dps{x_k\leftarrow \frac{1}{a_{kk}}\left[
                    a_{k,n+1}-\sum_{j=k+1}^na_{kj}x_j
                                               \right]}$\;
}}
\end{algorithm}



## Cas singulier et calcul du rang d'une matrice

\subsubsection{Pivotage total} Si, à l'itération $k$ de la méthode de Gauss
avec pivot total,
$$
\max\left\{\left|a^{(k)}_{ij}\right|\ ; \ i,j\in [\![k,n]\!]
    \right\}=0
$$
(donc tous les éléments de la sous-matrice sont nuls), on peut affirmer que
le rang de la matrice ${\bf A}$ est égal à $k-1$. S'il existe un élément 
$b^{(k)}_i$, $k\le i\le n$, du second membre différent de zéro, alors le
système n'a pas de solution.

\subsubsection{Pivotage partiel} Si
$$
\max\left\{\left|a^{(k)}_{ik}\right|\ ; \ \ i\in [\![k,n]\!]
    \right\}=0
$$
on peut seulement affirmer que la colonne $k$ est linéairement dépendante des
$k-1$ premières. Cela implique que $\mathrm{rang}({\bf A})<n$ et que le système n'a
probablement pas de solution. Toutefois, on ne peut en être sûr que sur
le test du pivot total nul.

\begin{rem}
La méthode de Gauss et le test du pivot total nul apparaissent donc comme 
la meilleure stratégie pour calculer le rang d'une matrice.
\end{rem}


# Facteurs LU d'une matrice non singulière

Quand on a plusieurs systèmes linéaires à résoudre avec la même matrice 
et des seconds membres différents, on a intérêt lors de la première 
résolution à garder les coefficients des pivotages successifs en mémoire.
Cela correspond à garder la factorisation $LU$\index{factorisation!LU}
de la matrice. En effet,
chaque pivotage peut être représenté par une matrice élémentaire qui
ne diffère de l'identité que par une sous-colonne.

Reprenons le pivotage de la matrice ${\bf A^{(k)}}$ à l'étape $k$. Soit ${\bm \eta_k}$
le vecteur de $\mathbb R^{n-k}$ dont les composantes sont
$\eta_{ik}=a^{(k)}_{ik}/a^{(k)}_{kk}$. On a donc (en supposant 
$a^{(k)}_{kk}\ne 0$)
$$
{\bf A^{(k+1)}}={\bf E_kA^{(k)}}
$$
où ${\bf E_k}$ est la matrice élémentaire suivante
\vskip 5pt
$$
  \left(
    \begin{array}{r@{}c|c@{}c}
   &  \begin{smallmatrix}
      \bovermat{k-1}{1 & & 0}\\
       &\ddots&\\
      0 & & 1\rule[-1ex]{0pt}{2ex}
     \end{smallmatrix} & \mbox{\huge0} & \rlap{\kern5mm k}\\\hline
   &  \mbox{\huge0} & 
     \begin{smallmatrix}\rule{0pt}{2ex}
     \vertbar&&&0\\
      -{\Large \bm \eta_k}& &  \\
       \vertbar&\ddots&\\
     & & & 1
     \end{smallmatrix}  & \rlap{\kern5mm n-k}
    \end{array} 
    \right)
$$


On a donc, après $n-1$ pivotages
\begin{equation}
{\bf A^{(n-1)}}={\bf E_{n-1}E_{n-2}}\cdots {\bf E_1 A}.\label{MatriceU}
\end{equation}
Notons ${\bf U}$ la matrice triangulaire supérieure ${\bf A^{(n-1)}}$ et
réécrivons la relation (\ref{MatriceU})
$$
{\bf A}={\bf E_1^{-1}}\cdots {\bf E_{n-1}^{-1} U}
$$ 
On vérifie que la matrice inverse ${\bf E_k^{-1}}$ a la même forme que
${\bf E_k}$ avec les éléments de la sous-colonne $k$ changés de signe et que les 
produits ${\bf E_k^{-1}E_{k+1}^{-1}}$ s'effectuent sans calcul en accolant les
vecteurs ${\bm\eta_k}$ et ${\bm \eta_{k+1}}$ dans les colonnes $k$ et $k+1$. On peut donc
écrire ${\bf A}={\bf LU}$ où ${\bf L}$ est une matrice triangulaire inférieure dont les éléments
diagonaux sont égaux à 1 et les éléments sous la diagonale sont
$$
l_{ij}=\eta_{ij}.
$$

\exemple{
Soit ${\bf A} = \begin{pmatrix}1&4&7\\2&5&8\\3&6&10\end{pmatrix}$. 
\begin{enumerate}
\item On pivote tout d'abord selon $a_{11}$. La matrice de transformation élémentaire est  \mbox{${\bf E_1} = \begin{pmatrix}1&0&0\\\textcolor{red}{-2}&1&0\\\textcolor{red}{-3}&0&1\end{pmatrix}$}, avec ${\bm \eta_1} = \begin{pmatrix} \textcolor{red}{2}\\\textcolor{red}{3}\end{pmatrix}$ et $ {\bf E_1A} =\begin{pmatrix}1&4&7\\0&-3&-6\\0&-6&-11\end{pmatrix}$
\item On pivote ensuite selon $a^{(2)}_{22}$. La matrice de transformation élémentaire est  \mbox{${\bf E_2} = \begin{pmatrix}1&0&0\\0&1&0\\0&\textcolor{blue}{-2}&1\end{pmatrix}$}, avec ${\bm \eta_2} = \begin{pmatrix} \textcolor{blue}{2}\end{pmatrix}$ et $ {\bf E_2E_1A} =\begin{pmatrix}1&4&7\\0&-3&-6\\0&0&1\end{pmatrix}={\bf U}$
\item Finalement ${\bf L} = \begin{pmatrix}1&0&0\\\textcolor{red}{2}&1&0\\\textcolor{red}{3}&\textcolor{blue}{2}&1\end{pmatrix}$
\end{enumerate}
}

Observons que les éléments non diagonaux de ${\bf L}$ peuvent être rangés
directement à la place des éléments de ${\bf A}$ correspondants. La matrice ${\bf A}$
est donc recouverte par sa factorisation ${\bf LU}$ et le coût de stockage est
en $n^2$.

Gr\^ace à cette factorisation (qui ne coûte donc pas plus cher que la
triangularisation), tout nouveau système linéaire ${\bf Ax}={\bf b'}$ peut être
résolu par la résolution de deux systèmes triangulaires (donc en $O(n^2)$ flops).
En effet, pour résoudre 
$$
{\bf LUx}={\bf b'}
$$
on résout d'abord 
$$
{\bf Ly}={\bf b'}
$$
puis 
$$
{\bf Ux}={\bf y}.
$$

\exemple{
Soit ${\bf A} = \begin{pmatrix}1&4&7\\2&5&8\\3&6&10\end{pmatrix}$. En appliquant l'algorithme de factorisation ${\bf LU}$ (cf. ci-dessus), on obtient 
$${\bf L} = \begin{pmatrix}1&0&0\\2&1&0\\3&2&1\end{pmatrix}\quad {\bf U} = \begin{pmatrix}1&4&7\\0&-3&-6\\0&0&1\end{pmatrix}$$
Si ${\bf b} = \begin{pmatrix}1\\1\\1\end{pmatrix}$, ${\bf Ly}={\bf b}$ donne ${\bf y} = \begin{pmatrix}1\\-1\\0\end{pmatrix}$ et ${\bf Ux}={\bf y}$ donne ${\bf x} =\frac{1}{3} \begin{pmatrix}-1\\1\\0\end{pmatrix}$
}

Quand aucun pivot nul n'est rencontré, ${\bf A}$ peut se mettre sous la forme ${\bf LU}$
et cette factorisation est unique. En effet, s'il existe deux factorisations
${\bf L_1U_1}$ et ${\bf L_2U_2}$ de ${\bf A}$, on a alors
${\bf L_1U_1}={\bf L_2U_2}$. Ce qui implique que ${\bf L_2^{-1}L_1}={\bf U_2U_1^{-1}}$ et le produit
de deux matrices inférieures (resp. supérieures) étant une matrice triangulaire
inférieure (resp. supérieure), ces produits sont nécessairement une matrice 
diagonale. C'est l'identité car $(l_1)_{ii}=(l_2)_{ii}=1$ pour tout $i$.

Dans le cas d'une stratégie de pivot partiel, si ${\bf P_k}$ est la matrice de 
permutation des lignes à l'itération $k$, on peut écrire
$$
{\bf A^{(k+1)}}={\bf E_kP_kA^{(k)}}.
$$
En fait, les différentes permutations peuvent être résumées dans la matrice
$$
{\bf P}={\bf P_{n-1}P_{n-2}}\cdots {\bf P_1}
$$
et on obtient la décomposition générale suivante :

\begin{theo} {Factorisation {\bf PA} = {\bf LU}}{}
Pour toute matrice ${\bf A}$ non singulière de taille $n$, il existe
une matrice de permutation ${\bf P}$, une matrice triangulaire inférieure ${\bf L}$ telle
que $l_{ii}=1$, pour tout $i$, et une matrice triangulaire supérieure ${\bf U}$, telles
que
$$
{\bf PA}={\bf LU}.
$$ 
\end{theo}

\medskip
L'algorithme~\ref{FactLU} montre les différentes étapes de la factorisation ${\bf LU}$
avec recherche du pivot partiel. En sortie, ${\bf A}$ contient les facteurs ${\bf L}$ et ${\bf U}$ de la matrice
et $\sigma$ les permutations de lignes éventuelles. 

\begin{algorithm}[H]

\Entree{La matrice ${\bf A}$}
\Sortie{Les facteurs ${\bf L}$ et ${\bf U}$, les permutations $\sigma$}
\Deb{
$\sigma(i)=i$, $i\in[\![1,n]\!]$ (initialisation du vecteur des permutations)\;
\Pour{$k=1,\ldots, n-1$}{
         \textit{Recherche du pivot}\;
         $c_p\leftarrow |a_{kk}|,\: i_p\leftarrow k$\;
         \Pour{$i=k+1,\ldots,n$}{
         \Si{$ |a_{ik}|> c_p$}{
           $c_p\leftarrow |a_{ik}|,\: i_p\leftarrow i$
         }
         }
         \textit{Permutation}\;
         \Si{$i_p\ne k$}{
     
         $\sigma(k)=i_p$, $\sigma(i_p)=k$
         }
         \textit{Pivotage}\;
    \Pour{$i=k+1,\ldots, n$}{
        \textit{Remplissage de la colonne} $k$ 
              \textit{par les coefficients} $\eta_{ik}$ \;
        $\dps{a_{ik}\leftarrow \frac{a_{ik}}{a_{kk}}}$\;
        \textit{Modification des lignes qui n'ont pas encore été ligne-pivot}\;
    \Pour{$j=k+1,\ldots,n+1$}{
         $\dps{a_{ij} \leftarrow  a_{ij}- a_{ik}a_{kj}}$\;
        }
    }
}}
\caption{Factorisation $LU$}
\label{FactLU}
\end{algorithm}



Si $\sigma_i=j$ alors les lignes
$i$ et $j$ ont été permutées. Les permutations doivent être repercutées sur le second
membre lors de la résolution de ${\bf Ly}={\bf b}$ (algorithme~\ref{SDI}).

\begin{algorithm}
\Deb{
\textit{1. Substitutions directes} $Ly=b$\;
$\dps{x_1\leftarrow {b_{\sigma_1}}}$\;
\Pour{$k=2,\ldots,n$}{
$\dps{x_k\leftarrow b_{\sigma_k}-\sum_{j=1}^{k-1}a_{kj}x_j}$\;
}
\textit{2. Substitutions inverses} $Ux=y$\;
$\dps{x_n\leftarrow \frac{x_{n}}{a_{nn}}}$\;
\Pour{$k=n-1,\ldots,1$}{
    $\dps{x_k\leftarrow \frac{1}{a_{kk}}\left[
                    x_{k}-\sum_{j=k+1}^na_{kj}x_j
                                               \right]}$\;
}}
\caption{Substitutions directes/inverses}
\label{SDI}
\end{algorithm}







## Cas particuliers


\subsubsection{Matrices symétriques} Dans ce cas, ${\bf U}$ peut s'écrire ${\bf U}={\bf DL^\top} $ où 
${\bf D}$ est la matrice diagonale contenant les pivots successifs.  On a donc
la factorisation ${\bf A}={\bf LDL^\top }$.\index{factorisation!${\bf LDL^\top} $} La complexité de 
l'algorithme est alors de $n^3/6$ flops (cf. exercice \ref{exo29}).

\subsubsection{Matrices bandes} Ce sont des matrices symétriques telles que 
$a_{ij}=0$ pour $|i-j|>p$ ($p$ est la largeur de bande de la matrice, $p<n$). 
Ces matrices interviennent couramment dans la discrétisation d'équations 
différentielles.  Il est alors facile de montrer que les facteurs ${\bf LU}$ respectent
la bande. On a alors intérêt de stocker la matrice (et ses facteurs ${\bf LU}$) sous
la forme de tableau à $n$ lignes et $p$ colonnes et la complexité est en $np^2/2$
flops.

\subsubsection{Matrices symétriques définies positives} Les matrices définies positives possèdent une factorisation unique ${\bf LDL^\top} $ avec des pivots
successifs strictement positifs. La factorisation peut s'effectuer directement sans
pivotage par identification terme à terme en $n^3/6$ flops par l'algorithme de Cholesky. 


# Autres applications

## Calcul de l'inverse d'une matrice
On a vu précédemment que la résolution d'un système linéaire ne nécessite pas le
calcul explicite de l'inverse d'une matrice. Quand on a besoin néanmoins de la 
calculer, on peut procéder de la manière suivante, basée sur la factorisation ${\bf LU}$
de la matrice:
\begin{itemize}
\item Calculer les facteurs ${\bf LU}$ de la matrice: ${\bf PA}={\bf LU}$
\item Résoudre les $n$ systèmes linéaires ${\bf LUx^i}={\bf Pe_i}$, où ${\bf e_i}$, $i\in [\![1, n]\!]$,
est le $i$-ème vecteur de la base canonique de $\mathbb R^n$. La solution ${\bf x^i}$ est la
$i$-ème colonne de ${\bf A^{-1}}$. 
\end{itemize}

Le coût total apparent est de $n^3/3+n^3=4n^3/3$ flops.
Mais on peut montrer que, gr\^ace à la structure particulière des seconds membres des
systèmes linéaires successifs, le coût réel n'est que de $n^3$ flops.

Une approche équivalente couramment utilisée, mais qui ne passe pas par le calcul 
des facteurs ${\bf LU}$, est la méthode dite de Gauss-Jordan qui consiste à pivoter
complètement le système paramétré
$$
{\bf Ax-y}={\bf 0}.
$$
On pivote cette fois sur la colonne entière de façon à transformer le système
en un système diagonal ${\bf x}-{\bf A^{-1}y}={\bf 0}$. On peut observer que cette technique consiste
à effectuer en parallèle à partir de la matrice identité les pivotages nécessaires
à la transformation de ${\bf A}$ en la matrice identité.


## Calcul du déterminant


Le calcul du déterminant d'une  matrice $2\times 2$ est bien connu:
\begin{equation}
\det\left[\begin{array}{rr} a & b \\ c & d\end{array}\right]=ad-bc \label{DetMat2x2}
\end{equation}
ce qui permet de déterminer explicitement l'inverse d'une matrice $2\times 2$:
\begin{equation}
\left[\begin{array}{rr} a & b \\ c & d\end{array}\right]^{-1}=\frac{1}{ad-bc}
\left[\begin{array}{rr} d & -b \\ -c & a\end{array}\right]. \label{InvMat2x2}
\end{equation}

La généralisation des formules (\ref{DetMat2x2})-(\ref{InvMat2x2}) au cas des matrices
$n\times n$ conduit aux fameuses formules de Cramer que nous ne reproduirons pas
ici car elles ont une complexité exponentielle, ce qui les rend impraticables. A titre d'exemple, pour calculer le déterminant d'une matrice
$20\times 20$ par la formule de Cramer il faut à peu près 15400 ans de calcul sur une
machine de 100 Mips (soit $10^8$ instructions par seconde). Avec la méthode des pivots
le coût n'est que de $3\cdot 10^{-5}$ secondes! 

En pratique on calculera le  déterminant après pivotage:
$$
\det({\bf A})=(-1)^p\prod_{i=1}^nu_{ii}
$$
où les $u_{ii}$ ($1\le i\le n$) sont les pivots et $p$ le nombre de permutations
effectuées au cours de la factorisation.






