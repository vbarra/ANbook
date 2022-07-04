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
# Projections orthogonales

On s'intéresse ici aux mécanismes de calcul de la projection orthogonale d'un vecteur de \mathbb R$^n$ sur un sous-espace. Les résultats fondamentaux d'existence, d'unicité et de caractérisation de la projection orthogonale d'un point sur un sous-espace, valables dans des espaces plus généraux que \R$^n$, sont rappelés ci-après. Dans tout le chapitre,
la norme vectorielle utilisée, et notée $\|\cdot\|$, est la norme euclidienne.

````{prf:theorem} Théorème de projection
Soit $L$ un sous-espace vectoriel de \mathbb R$^n$. Étant donné un vecteur ${\bf y}\in\mathbb R^n$, il existe un unique vecteur ${\bf p}\in L$, appelé projection orthogonale
de ${\bf y}$ sur $L$, tel que :

$(\forall {\bf x}\in L)\ \|{\bf y}-{\bf p}\|\leq \\{\bf y}-{\bf x}\|$
Une condition nécessaire et suffisante pour que ${\bf p}\in L$ soit la projection orthogonale de ${\bf y}$ sur $L$ est ${\bf y}-{\bf p}\in L^\bot$
````
```{index} Projection orthogonale
```


```{margin} 
![](./images/projection.png)
```

\textsc{Démonstration.} \\
\textit{1) Existence:} On peut supposer que ${\bf y}\notin L$ et on introduit 
$$\varepsilon=\inf_{{\bf x}\in L}\|{\bf y}-{\bf x}\|>0.$$ On veut trouver ${\bf p}\in L$ tel que $\|{\bf p}-{\bf y}\|=\varepsilon$. Soit une suite
$\{{\bf p_i}\}$ de vecteurs de $L$ tels que $\|{\bf y}-{\bf p_i}\|\rightarrow\varepsilon$. On a alors (loi du parallélogramme)
$$
\|({\bf p_j}-{\bf y})+({\bf y}-{\bf p_i})\|^2+\|({\bf p_j}-{\bf y})-({\bf y}-{\bf p_i})\|^2=2\|{\bf p_j}-{\bf y}\|^2+2\|{\bf y}-{\bf p_i}\|^2.
$$
Soit
$$
\|{\bf p_j}-{\bf p_i}\|^2=2\|{\bf p_j}-{\bf y}\|^2+2\|{\bf y}-{\bf p_i}\|^2-4\|{\bf y}-({\bf p_i}+{\bf p_j})/2\|^2.
$$
Comme $({\bf p_i}+{\bf p_j})/2\in L$, on a 
$$
\|{\bf y}-({\bf p_i}+{\bf p_j})/2\|^2\ge \varepsilon^2,
$$
ce qui permet de majorer
$$
\|{\bf p_j}-{\bf p_i}\|^2\le 2\|{\bf p_j}-{\bf y}\|^2+2\|{\bf y}-{\bf p_i}\|^2-4\varepsilon^2.
$$
Par passage à la limite, on déduit que $\|{\bf p_j}-{\bf p_i}\|\rightarrow 0$. Ce qui prouve que la suite des 
${\bf p_i}$ est une suite de Cauchy\footnote{Une suite de Cauchy est une suite de vecteurs $\{{\bf x_k}\}$
qui satisfait $\lim_{k,\ell\rightarrow+\infty}\|{\bf x_k}-{\bf x_\ell}\|=0$. Toute suite convergente est bien sûr de Cauchy
mais l'inverse n'est pas toujours vrai dans les espaces vectoriels normés généraux. C'est toujours vrai dans 
les espaces de Hilbert (donc dans $\mathbb R^n$).}


% On prend un $l \in L$ et on applique 
% le théorème de Weierstrass\footnote{Si $K$ est un compact non vide et $f$ est
% continue, alors le problème $\min_{x\in K} f(x)$ admet une solution} \index{théorème!de Weierstrass} 
% sur le compact  $L \cap B(y, \|y-l\|)$ pour la fonction continue et bornée inférieurement $\|x-y\|$.\\


\textit{2) Orthogonalité:} Supposons qu'il existe ${\bf x}\in L$ tel que ${\bf x^T} ({\bf y}-{\bf p})=\delta \neq 0$. On peut supposer que $\|{\bf x}\|=1$. On prend ${\bf p_1}={\bf p}+\delta {\bf x}$ et on écrit : \[  \|{\bf y}-{\bf p_1}\|^2=\|{\bf y}-{\bf p}\|^2-2 ({\bf y}-{\bf p})^T (\delta {\bf x}) + \delta^2 = \|{\bf y}-{\bf p}\|^2- \delta^2 < \|{\bf y}-{\bf p}\|^2 \]
d'où la contradiction.\\
\textit{3) Unicité:} Supposons qu'il existe un ${\bf x}\ne {\bf p}$ qui réalise le minimum.
On applique le théorème de Pythagore :
\[ (\forall {\bf x} \in L)\ \|{\bf y}-{\bf x}\|^2=\|{\bf y}-{\bf p}+{\bf p}-{\bf x}\|^2  = \|{\bf y}-{\bf p}\|^2 + \|{\bf p}-{\bf x}\|^2  > \| {\bf y}-{\bf p}\|^2 \]
d'où la contradiction. 

%--------------------------------------------------------------------------
\subsection{Projection sur une droite passant par l'origine}
\label{drt}
Une droite qui passe par l'origine est un sous-espace de dimension 1. Soit ${\bf y}\in \mathbb R^n$ et $D$ le sous-espace de vecteurs engendrés par un vecteur ${\bf v}$ non nul : $$D=\{{\bf z}\in \mathbb R^n /{\bf z}=\lambda {\bf v},\lambda\in \mathbb R\}$$
Si ${\bf p}\in D$ est la projection orthogonale de ${\bf y}$ sur la droite, on peut écrire  : ${\bf y}={\bf p}+{\bf u}$, avec ${\bf p}=\lambda {\bf v}$ et ${\bf u^T}{\bf v}=0$. (figure \ref{F:projD})\\

On remarque que la dernière relation signifie que ${\bf u}\in D^\bot$ et donc que ${\bf u}$ est la projection orthogonale de ${\bf y}$ sur $D^\bot$.
Alors :
$${\bf v^T}  {\bf y}=\lambda {\bf v^T}  {\bf v} \Rightarrow \lambda=\frac{1}{{\bf v^T}{\bf v} }. {\bf v^T}  {\bf y}$$
soit
$${\bf p}=\frac{ {\bf v^T}  {\bf y}}{ {\bf v^T}  {\bf v}}.{\bf v}$$

Ainsi, cette expression qui exprime le fait que ${\bf p}$ a la même direction que ${\bf v}$ peut s'écrire différemment pour représenter la transformation qui transforme ${\bf y}$ en ${\bf p}$.\\
On vérifie que $({\bf v^T}  {\bf y}){\bf v} = ({\bf v} {\bf v^T}){\bf y}$ où l'on observe que $ {\bf v^T}  {\bf y}\in \mathbb R$ alors que ${\bf v} {\bf v^T}\in\mathcal{M}_n(\mathbb R)$  de rang 1 (toutes les colonnes sont des multiples de ${\bf v}$) qui { projette} l'espace $\mathbb R^n$ sur la droite $D$. La projection orthogonale sur une droite qui passe par l'origine est donc une transformation linéaire, de matrice :\index{matrice!de projection}
$${\bf P}=\frac{1}{\norme{{\bf v}}^2}{\bf v}{\bf v^T} $$
\begin{figure}[hbtp!]
\begin{center}
\begin{tikzpicture}[line cap=round,line join=round,>=triangle 45,x=1.0cm,y=1.0cm,scale=0.65]
\clip(-3.52,-2.76) rectangle (7.36,5.76);
\draw [->] (0.,0.) -- (0.,4.);
\draw [->] (0.,0.) -- (-2.84,-2.02);
\draw[-](0,0) -- (7.68,-3.28);
\draw [->,line width=1.6pt,blue] (0.,0.) -- (3.84,3.62);
\draw [->,line width=1.2pt,cyan] (0.,0.) -- (3.84,-1.64);
\draw [->,line width=1.6pt,red] (0.,0.) -- (1.536,-0.656);
\draw [->,line width=1.6pt,orange] (0.,0.) -- (1.54,3.62);
\draw [dash pattern=on 4pt off 4pt] (3.84,3.62)-- (3.84,-1.64);
\draw [dash pattern=on 4pt off 4pt] (3.84,3.62)-- (1.54,3.62);
\draw [->] (0.,0.) -- (5.34,-1.3);
\draw (-0.54,4.74) node[anchor=north west] { $x_3$ };
\draw (-3.2,-1.4) node[anchor=north west] { $x_1$ };
\draw (5.48,-0.96) node[anchor=north west] { $x_2$ };
\draw (3.84,3.62) node[anchor=north west] { $\textcolor{blue}{{\bf y}}$ };
\draw (3.84,-1.64) node[anchor=north west] { $\textcolor{cyan}{{\bf p}}$ };
\draw (1.536,-0.656) node[anchor=north west] { $\textcolor{red}{{\bf v}}$ };
\draw (1.54,3.62) node[anchor=north west] { $\textcolor{red}{{\bf u}}$ };
\draw (5.76,-2.1) node[anchor=north west] { $D$ };
\end{tikzpicture}

\end{center}
\caption{Projection sur une droite : exemple dans $\mathbb R^3$}
\label{F:projD}
\end{figure}

{\rem :
\begin{enumerate}
  \item Comme ${\bf p}={\bf Py}$, la projection orthogonale sur le sous-espace orthogonal $D^\bot$ est \\${\bf u}={\bf y}-{\bf p}=(\mathbb I-{\bf P}){\bf y}$, donc $\mathbb I-{\bf P}$ est la matrice de projection orthogonale sur $D^\bot$.
  \item Soit $\theta$ l'angle entre les directions des vecteurs ${\bf v}$ et ${\bf y}$. On a alors :
  $$cos(\theta)=\frac{{\bf v^T}{\bf  y}}{\norme{{\bf v}}\norme{{\bf y}}}$$
  et on en déduit l'inégalité de Schwarz : 
  $$(\forall {\bf x},{\bf y}\in \mathbb R^n)\quad  |{\bf x^T}  {\bf y}|\leq \norme{{\bf x}}\norme{{\bf y}}$$
\end{enumerate}
}
\exemple{
Dans $\mathbb R^2$, la projection sur l'axe des abscisses, dirigé par ${\bf v}=\begin{pmatrix}1 \\0\\\end{pmatrix}$, est effectuée par \mbox{${\bf P}= \begin{pmatrix}1&0\\0&0\end{pmatrix}$} et on a bien ${\bf p}={\bf Py}=\begin{pmatrix}y_1 \\0\\\end{pmatrix}$}

\subsection{Projection sur une droite ne passant pas par l'origine}
On utilise la représentation d'une droite comme un sous-espace affine parallèle à un sous-espace : 
$$D=\{{\bf z}\in \mathbb R^n / {\bf z}={\bf z_0}+\lambda {\bf v},\lambda \in \mathbb R\}$$
Soit ${\bf y}\in \mathbb R^n$ et ${\bf p}$ sa projection orthogonale sur la droite $D$. Alors :\\
${\bf y}={\bf p}+{\bf u}$ avec ${\bf p}={\bf z_0}+\lambda {\bf v}$ et ${\bf u^T}  {\bf v}=0$, d'où ${\bf p}={\bf z_0}+{\bf P}({\bf y}-{\bf z_0})=(\mathbb I-{\bf P}){\bf z_0}+{\bf P}{\bf y}$,\\
 où ${\bf P}$ est la matrice de projection sur le sous-espace engendré par ${\bf v}$ obtenu dans le paragraphe \ref{drt}.
 
\subsection{Projection sur un sous-espace}\label{ssev}
Soit ${\bf A}\in\mathcal{M}_{n,r}(\mathbb R)$  de rang $r$ (donc $r\leq n$). Considérons la projection ${\bf p}$ d'un vecteur  ${\bf y}\in\mathbb R^n$ sur le sous-espace image de ${\bf A}$ (figure \ref{F:projImA}).
\begin{figure}[hbtp!]
\begin{center}
\begin{tikzpicture}[line cap=round,line join=round,>=triangle 45,x=1.0cm,y=1.0cm,scale=0.5]
\clip(-3.52,-2.76) rectangle (7.36,5.76);
\filldraw[fill=blue!20!] (-3.5,-2) -- (3.5,-2) --  (5.5,1) -- (-0.5,1) -- (-3.5,-2);
\draw[-](0,0) -- (0,4);

\draw [->,line width=1.6pt,blue] (0.,0.) -- (3.84,3.62);
\draw [->,line width=1.2pt,cyan] (0.,0.) -- (3.84,-1.64);
\draw [->,line width=1.6pt,orange] (0.,0.) -- (0,3.62);
\draw [dash pattern=on 4pt off 4pt] (3.84,3.62)-- (3.84,-1.64);
\draw [dash pattern=on 4pt off 4pt] (3.84,3.62)-- (0,3.62);
\draw (-0.54,4.74) node[anchor=north west] { $Ker({\bf A^T})$ };

\draw (3.84,3.62) node[anchor=north west] { $\textcolor{blue}{{\bf y}}$ };
\draw (3.84,-1.64) node[anchor=north west] { $\textcolor{cyan}{{\bf p}}$ };
\draw (1.54,3.62) node[anchor=north west] { $\textcolor{red}{{\bf u}}$ };
\draw (-3.5,-2) node[anchor=north west] { $ Im({\bf A})$ };
\end{tikzpicture}

\end{center}
\caption{Projection sur un sous-espace}
\label{F:projImA}
\end{figure}

Alors ${\bf y}={\bf p}+{\bf u}$ avec ${\bf p}={\bf Ax},{\bf x}\in \mathbb R^r$ et ${\bf A^T}{\bf u}=0$ car ${\bf u}\in (Im({\bf A}))^\perp = Ker({\bf A^T})$.\\
On obtient alors : ${\bf A^T} {\bf y}={\bf A^T}{\bf Ax}$, et comme ${\bf A}$ est de rang plein, la matrice ($r\times r$) ${\bf A^T A}\in\mathcal{M}_{r}(\mathbb R)$ est inversible et :\\
${\bf p}={\bf Py}$ avec ${\bf P}={\bf A(A^T A)^{-1}A^T }$.\\
De plus, on retrouve la matrice de projection orthogonale $\mathbb I -{\bf P}$ sur le noyau de ${\bf A^T }$.

{\rem Si $r=1$ on retrouve la formule trouvée dans le paragraphe \ref{drt}.}
%-----------------------------------
\subsection{Matrices de projection}
%-----------------------------------

Les matrices de projection\index{matrice!de projection} sont des matrices
qui possèdent les deux propriétés suivantes : 
\begin{itemize}
  \item ${\bf P}={\bf P^T} $ \rm [symétrie]
  \item ${\bf P}^2={\bf P}$\rm  [idempotence]
\end{itemize}
Les matrices de projection sont en général singulières, puisqu'elles ramènent l'espace à un sous-espace de dimension plus petite. De plus, elles contractent les normes : $\norme{{\bf Py}}\leq \norme{{\bf y}}$.