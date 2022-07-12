#!/usr/bin/env python
# coding: utf-8

# # Transformations orthogonales
# ## Matrices orthogonales
# ````{prf:definition} Matrice orthogonale
# Une matrice carrée ${\bf H}$ est dite {orthogonale}
# si et seulement si ${\bf H^T}{\bf H}={\bf H}{\bf H^T} = \mathbb I$
# ````
# ```{index} matrice;orthogonale
# ```
# ```{index} orthogonale;matrice
# ```
# 
# Une matrice orthogonale est donc une matrice carrée dont les colonnes sont orthonormées. Les matrices de rotation, de symétrie, de permutation et l'identité sont des exemples de matrices orthogonales.
# 
# Une matrice orthogonale ${\bf H}$ est naturellement inversible par définition, et l'inverse est ${\bf H^{-1}}={\bf H^T} $.
# 
# ````{prf:property} Propriété fondamentale}
# Les transformations orthogonales
# sont des isométries, les normes (euclidiennes), les produits scalaires et les angles sont conservés : 
# 
# $({\bf H} \textrm{ orthogonale}) \Leftrightarrow (\forall {\bf x}\in \mathbb R^n)\\|{\bf Hx}\|=\|{\bf x}\|$
# ````
# ```{index} transformation;orthogonale
# ```
# ```{index} orthogonale;transformation
# ```
# 
# En effet, $\|{\bf Hx}\|^2=({\bf Hx)^T} ({\bf Hx})={\bf x^T H^T Hx}={\bf x^T x}=\|{\bf x}\|^2$.
# 
# Cette propriété entraîne une stabilité numérique des méthodes utilisant ces transformations. On les utilise principalement pour :
# - orthonormaliser un système de générateurs,
# - résoudre un systèmes aux équations normales,
# - triangulariser un système mal conditionné,
# - calculer les valeurs propres d'une matrice.
# 
# 
# 
# Soit ${\bf Q}\in\mathcal{M}_n(\mathbb R)$ orthogonale, de colonnes 
# ${\bf q_1},{\bf q_2},\ldots,{\bf q_n}$. On a alors, pour tout ${\bf x}\in\mathbb R^n$, une représentation unique sur la base 
# orthonormée telle que
# 
# ${\bf x}=\displaystyle\sum_{i=1}^n ({\bf q_i^T} {\bf x}){\bf q_i}$
# 
#  $({\bf q_i^T} {\bf x}){\bf q_i}$ est la projection orthogonale de ${\bf x}$ sur l'axe ${\bf q_i}$. Cette représentation se
# généralise aisement à une base orthonormée quelconque ${\bf q_1},{\bf q_2},\ldots,{\bf q_r}$ d'un sous-espace de dimension $r$. 
# 
# ## Orthogonalisation de Gram-Schmidt}
# ```{index} orthogonalisation de Gram-Schmidt
# ```
# ```{index} Gram-Schmidt;orthogonalisation
# ```
# A partir d'une famille de $p$ vecteurs linéairement indépendants de $\mathbb R^n$, représentés par une matrice \mbox{${\bf A}\in\mathcal{M}_{n,p}(\mathbb R)$} de rang $p$, on peut construire une famille $\{{\bf q_1}\cdots {\bf q_p}\}$, base orthonormée de $Im({\bf A})$. C'est un outil fondamental pour la résolution de systèmes surdéterminés.
# 
# L'idée générale est donc de construire une base orthonormée du sous-espace image d'un ensemble de vecteurs. L'intérêt numérique est que cette construction équivaut à triangulariser la matrice formée par ces vecteurs. L'algorithme {prf:ref}`GS` et la figure suivante présentent le procédé d'orthonormalisation de Gram-Schmidt qui, s'il est simple à comprendre, ne présente qu'un intérêt académique puisqu'il est très coûteux. En pratique, on utilise la méthode de factorisation QR qui permet d'atteindre le même objectif grâce à des transformations orthogonales élémentaires (rotations de Givens ou transformations de  Householder) numériquement stables et seulement deux fois plus chères que la méthode de Gauss.
# 
# Le principe de Gram-Schmidt est de calculer, pour $j\in[\![2,p]\!]$,  chaque vecteur ${\bf q_j}$ en soustrayant à  ${\bf  A_{\bullet,j}}$ ses projections orthogonales sur les $j-1$
#  premiers vecteurs de la base orthonormée déjà calculés, puis en normant le résultat.
# 
# 
# ```{prf:algorithm} Procédé d'orthonormalisation de Gram-Schmidt - version de base
# :label: GS
# **Entrée : **  ${\bf A}\in\mathcal{M}_{n,p}(\mathbb R)$ de rang $p$
# 
# **Sortie : ** ${\bf Q_1}\in\mathcal{M}_{n,p}(\mathbb R)$ à colonnes orthonormées, ${\bf R_1}\in\mathcal{M}_{p}(\mathbb R)$ triangulaire supérieure
# 
# 1. $r_{11} = \\{\bf A_{\bullet,1}}\|$
# 2. ${\bf q_{{1}}} = \frac{{\bf A_{\bullet,1}}}{r_{11}}$
# 3. Pour $j$=2 à $p$
#     1. ${\bf p_j}={\bf A_{\bullet,j}}$
#     2. Pour $i$=1 à $j-1$
#         1. $r_{ij}={\bf A_{\bullet,j}^T}{\bf q_{{i}}}$
#         2. ${\bf p_j}={\bf p_j}- r_{ij}{\bf q_{{i}}}$
#     3. $r_{jj} = \|{\bf p_j}\|$
#     4. ${\bf q_{{j}}} = \frac{{\bf p_j}}{r_{jj}}$
# ```
# 
# 
# ![](./images/gs.png)
# 
# 
# 
# On remarque (voir boucles de l'algorithme) que la matrice ${\bf R_1}$ est triangulaire supérieure. 
# 
# Gram-Schmidt construit donc une matrice ${\bf Q_1}\in\mathcal{M}_{n,p}(\mathbb R)$ à colonnes orthormées et une matrice ${\bf R_1}\in\mathcal{M}_{p}(\mathbb R)$ telles que ${\bf Q_1^T}{\bf A} = {\bf R_1}$ soit ${\bf A}$=${\bf Q_1R_1}$.
# 
# ````{prf:example}
# Soit ${\bf A} = \begin{pmatrix} 1&1\\1&0\\1&1\end{pmatrix}$. 
# 
# Alors ${\bf A} ={\bf QR}$ avec ${\bf Q} = \begin{pmatrix} \textcolor{blue}{\frac{1}{\sqrt{3}}}&\textcolor{magenta}{\frac{1}{\sqrt{6}}}\\\textcolor{blue}{\frac{1}{\sqrt{3}}}&\textcolor{magenta}{-\frac{\sqrt 2}{\sqrt 3}}\\\textcolor{blue}{\frac{1}{\sqrt{3}}}&\textcolor{magenta}{\frac{1}{\sqrt{6}}}\end{pmatrix}$ et ${\bf R} = \begin{pmatrix} \textcolor{red}{\sqrt{3}}&\textcolor{orange}{\frac{2}{\sqrt{3}}}\\0&\textcolor{cyan}{\sqrt{\frac{2}{3}}}\end{pmatrix}$. En effet : 
# 
# 1. $\textcolor{red}{r_{11}} = \|{\bf A_{\bullet,1}}\| = \textcolor{red}{\sqrt{3}}$
# 2. $\textcolor{blue}{{\bf q_{1}}}= \frac{1}{\sqrt{3}}{\bf A_{\bullet,1}}= \textcolor{blue}{\frac{1}{\sqrt{3}}\begin{pmatrix}1\\1\\1\end{pmatrix}}$
# 3. $\textcolor{orange}{r_{12}}={\bf A_{\bullet,2}^Tq_1}=\begin{pmatrix}1&0&1\end{pmatrix}.\frac{1}{\sqrt{3}}\begin{pmatrix}1\\1\\1\end{pmatrix} =  \textcolor{orange}{\frac{2}{\sqrt{3}}}$
# 4. ${\bf p_2}={\bf A_{\bullet,2}}- r_{12}{\bf q_{{1}}}=\begin{pmatrix}1\\0\\1\end{pmatrix}-\frac{2}{\sqrt{3}}.\frac{1}{\sqrt{3}}\begin{pmatrix}1\\1\\1\end{pmatrix}=\begin{pmatrix}\frac{1}{3}\\-\frac{2}{3}\\\frac{1}{3}\end{pmatrix}$
# 5. $\textcolor{cyan}{r_{22}} = \|{\bf p_2}\|= \textcolor{cyan}{\sqrt{\frac{2}{3}}}$
# 6. $\textcolor{magenta}{{\bf q_{2}}} = \frac{{\bf p_2}}{r_{22}} = \textcolor{magenta}{\begin{pmatrix}\frac{1}{\sqrt{6}}\\-\frac{\sqrt 2}{\sqrt 3}\\\frac{1}{\sqrt{6}}\end{pmatrix}}$ 
# ````
# 
# 
# Il est possible de compléter ${\bf q_1}\cdots {\bf q_p}$ en une base orthonormée de $\mathbb R^n$, en continuant la procédure de Gram-Schmidt avec $n-p$ vecteurs arbitraires, mais tels que les $n$ colonnes formées avec les ${\bf A_{\bullet,j}}$ soient linéairement indépendantes. Soit ${\bf Q_2}$ la matrice des $n-p$ derniers vecteurs orthonormés. On a alors bien :
# 
# ${\bf A^T} {\bf Q_2}={\bf R_1^T} {\bf Q_1^T} {\bf Q_2}=0$ ce qui montre que :
# 
# ${\bf A}={\bf QR}=[{\bf Q_1}\ {\bf Q_2}]\cdot\left [\begin{array}{c}{\bf R_1} \\{\bf 0} \\\end{array}\right]={\bf Q_1R_1}$.\\
# 
# **Les colonnes de ${\bf Q_1}$ forment une base orthonormée de $Im({\bf A})$, et les colonnes de ${\bf Q_2}$ forment 
# une base orthonormée de $Ker({\bf A^T})$.** 
# 
# Appliqué au problème des moindres carrés, le système aux équations normales s'écrit donc 
# 
# $\begin{eqnarray*}
# {\bf A^TAx}&=& {\bf A^Tb}\\
# {\bf (Q_1R_1)^T(Q_1R_1)x}&=& {\bf (Q_1R_1)^Tb}\\
# {\bf R_1^TQ_1^TQ_1R_1x}&=& {\bf R_1^TQ_1^Tb}\\
# {\bf R_1^TR_1x}&=& {\bf R_1^TQ_1^Tb}\\
# {\bf R_1x}&=& {\bf Q_1^Tb}
# \end{eqnarray*}$
# 
# La dernière simplification étant possible car ${\bf R_1}$ est inversible (trangulaire supérieure, les éléments de la diagonale étant des normes, donc strictement positifs). La solution du 
# problème des moindres carrés est solution du système triangulaire 
# 
# ${\bf R_1x}= {\bf Q_1^Tb}$
# et le calcul de l'erreur donne
# 
# $\|{\bf e}\|^2=\|{\bf b}\|^2-\displaystyle\sum_{i=1}^p  ({\bf b^T} {\bf q_i})^2$
# 
# ## Transformations de Householder
# ```{index} transformation;de Householder
# ```
# On peut interpréter la méthode de Gram-Schmidt comme une méthode de triangularisation de la matrice ${\bf A}$, au même titre que la méthode de Gauss. Il est possible de réorganiser les calculs en construisant des transformations élémentaires orthogonales qui effectuent cette triangularisation colonne par colonne (ou élément par élément). Les symétries de Householder
# et les rotations de Givens sont des exemples simples et intéressants de telles transformations, car elles conduisent à des algorithmes numériquement plus stables que la méthode de Gram-Schmidt.
# 
# ```{margin} 
# ![](./images/householder.png)
# ```
# ````{prf:definition} Matrice de Householder}
# Une matrice de Householder\index{matrice!de Householder}
# est une matrice carrée ${\bf H}$ qui s'écrit ${\bf H}=\mathbb I-2{\bf P}$, où ${\bf P}$ est
# la matrice de projection\index{matrice!de projection}
# sur la droite engendrée par un vecteur ${\bf v}$ non nul.
# ````
# On vérifie  que ${\bf H}$ représente une symétrie
# par rapport au sous-espace ${\bf v^\perp}$.
# ```{index} symétrie
# ```
# 
# 
# Le théorème suivant montre qu'il est toujours possible de trouver une matrice de Householder\index{matrice!de Householder}
# permettant de transformer un vecteur quelconque en un vecteur colinéaire à un vecteur donné 
# \begin{theo}{}{}
# \label{T:householder}
# Soient ${\bf f}$ et ${\bf e}$ deux vecteurs non colinéaires de $\mathbb{R}^n$; avec $\|{\bf e}\|_2=1$. Il est alors possible de trouver ${\bf u}\in \mathbb{R}^n$ tel que 
# \begin{enumerate}
#   \item $\|{\bf u}\|_2=1$
#   \item ${\bf H(u)f}=\alpha {\bf e}$
# \end{enumerate}
# \end{theo}
# 
# \textsc{Démonstration:} 
# Remarquons tout d'abord que si ${\bf H(u)}$ est une matrice de Householder, alors ${\bf H(u)f}={\bf f-2u(u^T f)}$ et $\|{\bf H(u)f}\|_2=\|{\bf f}\|_2$.\\
# Posons alors $\mid\alpha\mid=\|{\bf f}\|_2$. On cherche alors ${\bf u}$ tel que ${\bf H(u)f}=\alpha {\bf e}$, soit 
# \begin{align*}
# {\bf f-2u(u^T f)}&=\alpha {\bf e}\\
# {\bf u}&=\frac{1}{2{\bf u^T f}}({\bf f-}\alpha {\bf e})
# \end{align*}
# Si $\beta={\bf u^T f}$, en multipliant à gauche par ${\bf f^T} $ :
# $$2\beta^2=\alpha^2-\alpha {\bf f^T e}$$
# et $\beta$ existe si $\alpha^2-\alpha {\bf f^T e}>0$. Or l'inégalité de Cauchy-Schwarz
# \index{inégalité!de Cauchy-Schwarz} nous donne
# $$\mid {\bf f^T e}\mid\leq\|{\bf f}\|_2\|{\bf e}\|_2=\|\alpha\|$$
# et l'inégalité est de plus stricte par hypothèse (${\bf f}$ et ${\bf e}$ non colinéaires). Ainsi :
# $${\bf u}=\frac{1}{2\beta}({\bf f}-\alpha {\bf e)}$$ répond à la question.
# 
#  
# \begin{rem}
# Si ${\bf f}$ et ${\bf  e}$ sont colinéaires, ${\bf H}=\mathbb I$ ou ${\bf H}=\mathbb I-2{\bf ee^T }$ répondent à la question.
# \end{rem}
# 
# L'algorithme d'orthonormalisation  de ${\bf A}$ par matrices de Householder  opère alors colonne par colonne, et transforme itérativement ${\bf A}$ en une matrice triangulaire supérieure.
# 
# 
# On l'illustre dans la suite (algorithme \ref{A:HS}) dans le cas où ${\bf A}\in\mathcal{M}_n(\mathbb R)$ est de rang plein.
# 
# \begin{algorithm}[H]
# \Entree{${\bf A}\in\mathcal{M}_n(\mathbb R)$}
# \Sortie{${\bf Q}\in\mathcal{M}_n(\mathbb R)$ orthogonale, ${\bf R}\in\mathcal{M}_n(\mathbb R)$ triangulaire supérieure}
# \Deb{
# ${\bf A^{(1)}}={\bf A}$\\
# \Pour {$j$=1 à $n-1$}
# {
# (i) Soit ${\bf f_j}\in \mathbb R^{n-j+1}$ le vecteur commençant à l'élément $(j,j)$ de ${\bf A^{(j)}}$\\
# (ii) On construit  (théorème \ref{T:householder}) ${\bf {\tilde H^{(j)}}}\in\mathcal{M}_{n-j+1}(\mathbb R)$ telle que ${\bf \tilde{H^{(j)}}f_j} = \norme{{\bf f_j}} {\bf e^{(j)}_1}$, ${\bf e^{(j)}_1}$ premier vecteur de la base canonique de $\mathbb R^{n-j+1}$\\
# (iii) On construit \[{\bf H^{(j)}} =
# \left (
# \begin{array}{lll|l}
# \bovermat{j-1}{1 &\cdots &0&}\\
#   \vdots&\ddots &\vdots &{\bf 0}\\
# 0 & \cdots& 1& \\
# \hline
#  & {\bf 0}& & {\bf \tilde{H^{(j)}}}\\
#  
# \end{array}
# \right )\in\mathcal{M}_n(\mathbb R)
# \] 
# 
# (iv) On calcule ${\bf A^{(j+1)}} = {{\bf H^{(j)}}\bf A^{(j)}}$
# }
# ${\bf R}={\bf A^{(n-1)}}$ et ${\bf Q} =  {\bf {H^{(1)}}^T}{\bf {H^{(2)}}^T} \cdots  {\bf {H^{(n-1)}}^T}$
# }
# \caption{Factorisation {\bf QR} par matrices de Householder}
# \label{A:HS}
# \end{algorithm}
# 
# A l'issue des $n-1$ itérations, on a effectué les produits ${\bf H^{(n-1)}} \cdots {\bf H^{(2)}} {\bf H^{(1)}}$ pour obtenir une matrice triangulaire supérieure ${\bf R}\in\mathcal{M}_n(\mathbb R)$ à partir de ${\bf A}$. Donc :
# $${\bf H^{(n-1)}} \cdots {\bf H^{(2)}} {\bf H^{(1)}}{\bf A} = {\bf R}$$
# De plus les ${\bf H^{(j)}}, j\in[\![1, n-1]\!]$ sont orthogonales donc
# $${\bf A} ={\bf {H^{(1)}}^T}{\bf {H^{(2)}}^T} \cdots  {\bf {H^{(n-1)}}^T} {\bf R}$$
# Le produit des matrices ${\bf{ H^{(j)}}^T}$ est également une matrice orthogonale, et on pose 
# $${\bf Q} = {\bf {H^{(1)}}^T}{\bf {H^{(2)}}^T} \cdots  {\bf {H^{(n-1)}}^T}$$ 
# pour finalement obtenir 
# $${\bf A} = {\bf QR}$$
# 
# On illustre cet algorithme sur les deux premières itérations :
# \begin{enumerate}
# \item $j$=1 : 
# \begin{maliste}
# \item On construit ${\bf H^{(1)}}$ telle que ${\bf H^{(1)}A^{(1)}_{\bullet,1}}$ =${\bf e_1}$, premier vecteur de la base canonique de $\mathbb R^n$
# \item ${\bf A^{(2)}}={\bf H^{(1)}A^{(1)}} = \begin{pmatrix}{\norme{{\bf A^{(1)}_{\bullet,1}}}}&a^{(2)}_{12}&*&\cdots &*\\0&a^{(2)}_{22}&*&\cdots &*\\0&\vdots&\vdots&\vdots&\vdots\\0&a^{(2)}_{n2}&*&\cdots &*\end{pmatrix}$
# \end{maliste}
# \item$ j$=2 : 
# \begin{maliste}
# \item soit ${\bf f_2} = \begin{pmatrix}a^{(2)}_{22}\\\vdots\\a^{(2)}_{n2}\end{pmatrix}\in \mathbb R^{n-1}$
# \item On construit ${\bf \tilde{H^{(2)}}}\in\mathcal{M}_{n-1}(\mathbb R)$ telle que ${\bf \tilde{H^{(2)}}f_2} = \norme{{\bf f_2}} {\bf e^{(2)}_1}$, ${\bf e^{(2)}_1}$ premier vecteur de la base canonique de $\mathbb R^{n-1}$
# \item On construit
# ${\bf H^{(2)}} =\left (
# \begin{array}{l|l}
# 1&{\bf 0}\\
# \hline
#  {\bf 0} & {\bf \tilde{H^{(2)}}}\\
# \end{array}
# \right )\in\mathcal{M}_n(\mathbb R)
# $ telle que $${\bf A^{(3)}}={\bf H^{(2)}A^{(2)}}\begin{pmatrix}{\norme{{\bf A^{(1)}_{\bullet,1}}}}&\ast&\ast&\cdots &\ast\\0&\norme{{\bf f_2}}&a^{(3)}_{23}&\cdots &\ast\\0&0&a^{(3)}_{33}&\vdots&\vdots\\0&0&a^{(3)}_{n3}&\cdots &*\end{pmatrix}$$
# \end{maliste}
# \end{enumerate}
# 
# \aretenir{
#   \begin{enumerate}
#   \item Projections : définitions, propriétés.
#   \item Moindres carrés : savoir résoudre un problème aux moindres carrés, en passant par le système aux équations normales
#   \item Savoir décomposer une matrice A sous la forme QR, par Gram Schmidt
#   \item Interpréter géométriquement cette décomposition.
# \end{enumerate}
# }
