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
# Elimination de Gauss
## Exemple introductif
\exemple{


\label{UnExemple}
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
les variables calculées dans les équations ayant servi à l'élimination :\\
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
}
On remarque qu'à une étape donnée, l'élimination d'une variable peut se faire
dans n'importe quelle équation, à condition, bien s\^ur, que cette équation
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
On remarque que le calcul de $x_k$ co\^ute $n-k$ flops et une division. Le
co\^ut total de l'algorithme est donc de $$1+2+\cdots+n-1=\frac{n(n-1)}{2},$$
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
Le co\^ut de l'algorithme de substitution (\ref{FwdSub1})-(\ref{FwdSub2}) est
le même que celui de substitution inverse (\ref{BckSub1})-(\ref{BckSub2}).

## Méthode de Gauss (ou du pivot) pour les systèmes linéaires


Montrons d'abord que la technique d'élimination décrite au paragraphe \ref{UnExemple}
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
l'étape $k$ co\^ute $(n-k)(n-k+1)$ flops et $n-k$ divisions. On en déduit le 
co\^ut total de l'algorithme de Gauss:
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
(la recherche du plus grand élément sur une liste de $p$ nombres co\^ute
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
probablement pas de solution. Toutefois, on ne peut en être s\^ur que sur
le test du pivot total nul.

\begin{rem}
La méthode de Gauss et le test du pivot total nul apparaissent donc comme 
la meilleure stratégie pour calculer le rang d'une matrice.
\end{rem}


