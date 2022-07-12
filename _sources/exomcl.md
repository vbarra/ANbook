\begin{exo}\rm
Trouver la meilleure approximation, $\bar x$, au sens des moindres carrés, du
système
$$
3x=10,\qquad 4x= 5.
$$
Déterminer le carré de l'erreur $e^2$ et montrer que le vecteur erreur
$(10-3\bar x\;  5-4\bar x)^T $ est orthogonal à $(3\; 4)^T $
\end{exo}


\begin{exo}\rm
Trouver la droite $f(t)=at+b$ qui approche le mieux (au sens des moindres
carrés) les mesures : $$ f(0)=0,\quad f(1)=1,\quad f(3)=2,\quad f(4)=5. $$
\end{exo}


\begin{exo}\rm
Soit $$P=\frac{1}{\parallel v\parallel^2}vv^T ,$$ la matrice de 
projection sur la droite de direction $v$.

Interpréter géométriquement la matrice $Q=I-2P$. Montrer que $Q^2=I$ et
$Q^T =Q$ (\textit{i.e.} $Q$ est orthogonale). Soit $y=(y_1, y_2)\in\mbb R^2$.
Déterminer $v$ pour que la seconde coordonnée de $z=Qy$ soit nulle.
\end{exo}

\begin{exo}\label{exo44}\rm
Projeter le vecteur $b=(0\; 3\; 0)^T $ sur les droites de directions 
respectives 
$$
a^1=\left(\begin{array}{r}2/3\\ 2/3\\ -1/3\end{array}\right),\quad
a^2=\left(\begin{array}{r}-1/3\\ 2/3\\ 2/3\end{array}\right).
$$

Trouver la projection de $b$ sur le plan engendré par $a^1$ et $a^2$.
\end{exo}

\begin{exo}[Examen novembre 2005]\rm
 On se donne
$$
A_1=\left[\begin{array}{cc}1 & -1\\ 1 & 0\\ 1 & 1\\ 1 & 2\end{array}
  \right]\quad\mbox{et }\ 
b=\left[\begin{array}{c}1 \\ 1\\ 0\\ 1 \end{array}
  \right]
$$

a) Mettre $A_1$ sous la forme $Q_1R_1$ où $Q_1$ est une matrice $4\times 2$ de colonnes
orthonormées $q_1$ et $q_2$, et $R_1$ une matrice $2\times 2$ triangulaire supérieure.

b) Calculer la solution aux moindres carrés du système incompatible
$$
A_1x_1=b
$$
et l'erreur $e_1$ correspondante.

c) On rajoute la colonne $a_3=[1\ 0\ 1\ 4]^T $ à la matrice $A_1$ pour avoir
$A_2=[A_1\ a_3]$. Mettre $A_2$ sous la forme $Q_2R_2$ et calculer la solution aux moindres
carrés du système
$$
A_2x_2=b
$$ 
Calculer l'erreur $e_2$ correspondante et montrer que $e_2=e_1^2-(q_3^T b)^2$,
où $q_3$ est la troisième colonne de $Q_2$. Justifier cette diminution dans le cas général.
\end{exo}

\begin{exo}[Examen février 2001]\rm
 On se propose de déterminer les coefficients du polynôme
$$
p(x)=\sum_{k=0}^na_kx^k
$$
qui approche <<au mieux>> une fonction $y=f(x)$, connue à travers $m+1$ 
couples $(x_i,y_i)$, $i=0,\ldots,m$; avec $m\ge n$. On supposera que tous les
points de mesure $x_i$, $i=0,\ldots,m$, sont {\it distincts}. 
Pour la suite on pose
${\bf x}=(x_0,\ldots,x_m)$ et ${\bf y}=(y_0,\ldots,y_m)$.

Pour ${\bf a}=(a_0,a_1,\ldots,a_n)\in \mbb R^{n+1}$, on définit
$$
J({\bf a})=\sum_{i=0}^m\Bigl[y_i-p(x_i)\Bigr]^2.
$$
\begin{enumerate}
\item Caractériser la solution ${\bf a}^\ast$ qui minimise $J({\bf a})$
sur $\mbb R^{n+1}$. Quelle est l'interprétation du polynôme $p(x)$
quand $m=n$. Quelle est alors la valeur de $J({\bf a}^\ast)$.

\item Montrer qu'il existe une matrice $U$, rectangulaire à $(m+1)$ lignes
et $(n+1)$ colonnes, telle que, si ${\bf a}^\ast$ est un point stationnaire de
$J$, alors
$$
U^{T}U{\bf a}^\ast=U^{T}{\bf y}.
$$
\item On suppose que $n=1$, i.e. $p(x)=a_0+a_1x$. On pose
\begin{eqnarray*}
&&\bar{{\bf x}}=\frac{1}{m+1}\sum_{i=0}^mx_i,\qquad
           \bar{{\bf y}}=\frac{1}{m+1}\sum_{i=0}^my_i,\\[3pt]
&&\sigma({\bf x}^2)=\frac{1}{m+1}\sum_{i=0}^m(x_i-\bar{{\bf x}})^2,\qquad
  \sigma({\bf x},{\bf y})=\frac{1}{m+1}\sum_{i=0}^m
         \left(x_i-\bar{{\bf x}}\right)\left(y_i-\bar{{\bf y}}\right).
\end{eqnarray*}
Montrer que
$$
a_1=\frac{\sigma({\bf x},{\bf y})}{\sigma({\bf x}^2)},\qquad
a_0=\bar{{\bf y}}-a_1\bar{{\bf x}} 
$$
et en déduire que le <<point moyen>> $(\bar{{\bf x}},\bar{{\bf y}})$ 
appartient à la droite $y=a_0+a_1x$.
\end{enumerate}
\end{exo}

\begin{exo}[Examen décembre 1996]\rm 
 On se donne 
$$
A = \left[\begin{array}{ccc}1 & 2 & -1 \\ 2 & 1 &  4 \\ 1 & 1 &  1\end{array}\right]
$$ 
et 
$$
b = \left[\begin{array}{c}1 \\ 1 \\1 \end{array}\right].
$$

\begin{enumerate}
\item Montrer que le système $Ax = b$ est incompatible.
\item Donner une CNS pour que $x \in \mbb R^3$
réalise le minimum de $E(x) = \| Ax - b \|_2^2$.
\item En déduire que la solution des moindres carrés n'est pas unique.  
\item Construire une base orthonormée de $ImA$ par le procédé
d'orthonormalisation de Gram-Schmidt. En déduire l'existence
d'une matrice $Q$, de format $3\times 2$, formée de vecteurs colonnes
orthonormés, et d'une matrice $R$ triangulaire supérieure de format
$2\times 3$, telles que $A = QR$.
\item Montrer que $RR^T $ est une matrice inversible. En déduire
une expression explicite pour l'ensemble des solutions des moindres
carrés de la forme $x = x_0 + \alpha u$, où $\alpha \in \mbb R$, et
tel que $x_0$ et $u$ soient des vecteurs orthogonaux.
% \item Exprimer la solution des moindres carrés de norme minimale en fonction
% de $Q$, $R$ et $b$. On ne demande pas de la calculer mais de raisonner avec les projecteurs.
\item Exprimer $\min\| Ax - b \|_2^2$ en fonction de $b$ et d'un 
vecteur orthogonal à $\mathrm{Im}Q$.  
\item Application numérique : calculer la solution des moindres carrés de norme minimale, ainsi que l'erreur résiduelle $\min\| Ax - b \|_2^2$. 
\end{enumerate} 
\end{exo}


\begin{exo}\rm
Etant donn\'ee une matrice $A (m\times n)$, on veut construire une matrice $M (m \times m)$ telle que 
\begin{itemize}
 \item $MA=S$, o\`u $S$ est triangulaire sup\'erieure $(m \times n)$
  \item $MM^T =\Delta^2$, o\`u $\Delta=\textrm{Diag}\{\delta_1,\ldots, \delta_m\}$ avec $\delta_i \neq 0, i=1,\ldots,m$
\end{itemize}

\begin{enumerate}
 \item Montrer que le calcul de $M$ permet d'obtenir la factorisation $QR$ de $A$ et que $Q=M^T \Delta^{-1}$, $R=\Delta^{-1}S$
 \item Analysons le cas $m=2$ : soient $x\in \mbb R^2$ et $\Delta={\textrm Diag}\{\delta_1, \delta_2\}$ ($\delta_i \neq 0$) donn\'es.
  \begin{enumerate}
   \item On d\'efinit
   \[ M_1 = \left[ \begin{array}{lr}
   \beta_1 & 1 \\ 1 & \alpha_1\end{array} \right] \]
   Supposer $x_2 \neq 0$ et calculer $M_1 x$ et $M_1\Delta^2 M_1^T $.
   
   Comment choisir $\alpha_1$ et $\beta_1$ de fa\c{c}on \`a ce que la deuxi\`eme composante de $M_1x$ soit nulle et que $M_1\Delta^2 M_1^T $
   soit diagonale ?
   
   Pour le choix pr\'ec\'edent, d\'eterminer $\gamma_1$ tel que 
   \[ M_1x = \left[ \begin{array}{c} (1+ \gamma_1)x_2 \\ 0 \end{array} \right] \ \ \mbox{et} \ \  
          M_1\Delta^2 M_1^T  = \left[ \begin{array}{lr} (1+ \gamma_1)\delta_2^2 & 0 \\ 
                                               0  & (1+ \gamma_1)\delta_1^2 
   \end{array} \right] \]
   
   \item Supposer $x_1 \neq 0$; on d\'efinit
   \[ M_2 = \left[ \begin{array}{lr}
   1 & \alpha_2 \\ \beta_2 & 1 \end{array} \right] \]
   Choisir $\alpha_2$ et $\beta_2$ de fa\c{c}on \`a ce que
    \[ M_2x = \left[ \begin{array}{c} (1+ \gamma_2)x_1 \\ 0 \end{array} \right] \ \ \mbox{et} \ \  
     M_2\Delta^2 M_2^T  = \left[ \begin{array}{lr} (1+ \gamma_2)\delta_1^2 & 0 \\ 0  & (1+ \gamma_2)\delta_2^2 
   \end{array} \right] \]
  et d\'eterminer $\gamma_2$.
  \end{enumerate}
 \item Soit maintenant $m$ un entier quelconque; d\'efinir les matrices $M_1(p,q)$ et $M_2(p,q)$ telles que
 \[ \left[ \begin{array}{lr} m_{pp} & m_{pq} \\ m_{qp} & m_{qq} \end{array} \right] = \left[ \begin{array}{lr} \beta_1 & 1 \\ 1 & \alpha_1  \end{array} \right]  \ \ \mbox{ou} \ \ 
 = \left[ \begin{array}{lr} 1 & \alpha_2 \\ \beta_2 & 1  \end{array} \right]  \]
 \[ \mbox{la composante q de} \ M_i(p,q)x\ \mbox{est nulle} \]
 
 \[ M_i\Delta^2 M_i^T  \ \mbox{diagonale} \]
 
%\newpage


 Les matrices $M_i$ sont appel\'ees {\sl matrices de Givens rapide}. On r\'esume ci-dessous l'algorithme de Givens rapide :
 
\begin{algorithm}
  $\delta_i=1,i=1,\ldots m$\;
  \Pour{$p=1,\ldots \min\{n,m-1\}$}{
   \Pour{$q=p+1,\ldots m$}{
   \Si{$a_{pq}\neq 0$} {
     \[ \alpha = - a_{pp}/a_{qp},\quad \beta = - \alpha \delta_q/\delta_p,\quad \gamma= - \alpha \beta \]
   }
   \Si{$\gamma \leq 1$} {
   \[ \left[ \begin{array}{lcr} a_{pp} & \cdots & a_{pn} \\
   a_{qp} & \cdots & a_{qn} \end{array} \right] \ \leftarrow  \left[ \begin{array}{lr} \beta &1 \\ 1 & \alpha \end{array} \right] \left[ \begin{array}{lcr} a_{pp} & \cdots & a_{pn} \\
   a_{qp} & \cdots & a_{qn} \end{array} \right] \]
   $\delta_p \ \leftarrow (1+ \gamma \delta_p)$\\
   $\delta_q \ \leftarrow (1+ \gamma \delta_q)$}
   \Sinon{
   \'echanger $\alpha$ et $\beta$\\
   \[ \alpha = 1/\alpha,\quad \beta= 1/\beta,\quad \gamma=1/\gamma \]\\
   \[ \left[ \begin{array}{lcr} a_{pp} & \cdots & a_{pn} \\
   a_{qp} & \cdots & a_{qn} \end{array} \right] \ \leftarrow  \left[ \begin{array}{lr} 1 & \alpha \\ \beta & 1 \end{array} \right] \left[ \begin{array}{lcr} a_{pp} & \cdots & a_{pn} \\
   a_{qp} & \cdots & a_{qn} \end{array} \right] \]\\
   $\delta_p \ \leftarrow (1+ \gamma \delta_p)$\\
    $\delta_q \ \leftarrow (1+ \gamma \delta_q)$\\
   }
   }
   }
\end{algorithm}

 \item Justifier l'algorithme et donner son co\^{u}t en nombre de flops. Comparer avec la m\'ethode de Householder vue en cours.
 
 \item Application num\'erique :
 r\'esoudre au sens des moindres carr\'es par la m\'ethode de Givens rapide le syst \`eme incompatible $Ax=b$ avec 
 \[ A= \left[ \begin{array}{lr} 1 & 4 \\ 2 & 5 \\ 3 & 6 \end{array} \right], b = \left[ \begin{array}{c} 7 \\ 8 \\ 9 \end{array} \right] \]
 
 \end{enumerate}
 
 \end{exo}