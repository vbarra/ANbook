#!/usr/bin/env python
# coding: utf-8

# # Moindres carrés linéaires
# 
# ## Modélisation
# On observe un système $\mathcal S$ qui produit, pour une sollicitation en entrée $t$ donnée, une sortie $Y$ (par exemple la vitesse d'un coureur en fonction du temps, la température d'une pièce mesurée chaque seconde,...). L'objectif est de construire un modèle théorique de $\mathcal S$, qui reproduise aussi fidèlement que possible le comportement du système.
# 
# On se donne donc un modèle (une fonction) $f_{\bf x}$, paramétrée par $n$ paramètres ${\bf x}=\left (x_1\cdots x_n\right )^T\in\mathbb R^n$, qui calcule pour chaque entrée $t$ la valeur $f_{\bf x}(t) = y$. 
# 
# Pour $m$ entrées $t_i, 1\leq i\leq m$, on souhaite donc que les valeurs $y_i=f_{\bf x}(t_i)$ calculées par le modèle soient les plus proches possibles des valeurs théoriques $Y_i$ mesurées sur $\mathcal S$. L'ajustement du modèle se fait à l'aide du vecteur de paramètres ${\bf x}$.
# 
# L'erreur entre le modèle théorique et la sortie mesurée pour l'entrée $t_i$ est $e_i=y_i-Y_i$ et la solution aux moindres carrés 
# consiste à choisir  ${\bf x}$ qui minimise la somme des carrés des erreurs. On pose donc le problème :
# 
# 
# Trouver ${\bf x}\in \mathbb R^n$ tel que $\displaystyle\sum_{i=1}^me_i^2$ soit minimale.
# 
# 
# Si la fonction $f_{\bf x}$ est  linéaire par rapport aux paramètres ${\bf x}$, alors le modèle s'écrit 
# $f_{\bf x}(t_i) = {\bf a_i}^T{\bf x}$ et on parle de problème au moindres carrés linéaires.
# 
# Soient ${\bf Y}\in \mathbb R^m$ le vecteur des valeurs $Y_i$ mesurées sur $\mathcal{S}$, ${\bf e}\in \mathbb R^m$ le vecteur dont les composantes sont les erreurs $e_i$ et  ${\bf A}\in\mathcal{M}_{mn}(\mathbb R)$ la matrice dont les lignes sont les ${\bf a_i^T} $. Le problème consiste alors à :
# 
# 
# Trouver ${\bf x}\in \mathbb R^n$ qui  minimise $\|{\bf Ax}-{\bf Y}\|^2$.
# 
# 
# Dans $\mathbb R^m$, il s'agit donc de trouver le point de l'image de ${\bf A}$ le plus proche au sens de la norme euclidienne du vecteur ${\bf Y}$ (qui a peu de chances d'appartenir à $Im({\bf A})$, car $m>>n$). L'unique solution (cf.  théorème \ref{T:proj}) est  la  projection orthogonale du vecteur ${\bf Y}$ sur le sous-espace $Im({\bf A})$. La solution a déjà été calculée au paragraphe \ref{ssev}, c'est la solution du système linéaire suivant, dit système aux équations normales : ${\bf A^T Ax}={\bf A^T Y}$
# 
# Si $rang({\bf A})=n$ (hypothèse raisonnable car $m>>n$ et les ${\bf a_i}$ dépendent des entrées $t_i$), ce système a une solution unique 
# 
# ${\bf x}={\bf (A^T A)^{-1}A^T Y}$
# 
# dite solution aux moindres carrés. La matrice ${\bf A^+}={\bf (A^T A)^{-1}A^T}$ (Attention !! ce n'est pas la matrice de projection) est la pseudo inverse de la matrice rectangulaire ${\bf A}$. Elle satisfait ${\bf A^+A}=\mathbb I$ et ${\bf AA^+}={\bf P}$.
# 
# 
# ````{prf:example} 
# Le système $\mathcal{S}$ produit pour les entrées $t_i =\{-1,0,1\}$ les sorties $Y_i=\{4,5,9\}$
# 
# On modélise $\mathcal{S}$ à l'aide d'un modèle affine $y=\alpha+\beta t$, où $\alpha$ et $\beta$ sont deux paramètres à déterminer pour minimiser au sens des moindres carrés l'erreur des trois mesures.
# 
# | Point | $t_i$ | $Y_i$ | $y_i$          |
# |-------|-------|-------|----------------|
# | $P_1$ | -1    | 4     | $\alpha-\beta$ |
# | $P_2$ | 0     | 5     | $\alpha$       |
# | $P_3$ | 1     | 9     | \alpha+\beta$  |
# 
# 
# Le système 
# \[{\bf Ax} = {\bf Y} :
# \left \{
# \begin{array}{ccccc}
# \alpha & - & \beta& = & 4 \\
# \alpha &  & & = & 5 \\
# \alpha & + & \beta& = & 9 \\
# \end{array}
# \right.
# \] 
# n'a bien sûr pas de solution. On construit les équations normales en multipliant 
# le système par 
# $${\bf A^T} =\left [\begin{array}{ccc}1 & 1 & 1\\-1 & 0 & 1\\\end{array}\right]$$
# Comme
# $
# {\bf A^T A}=\left[\begin{array}{cc}3 & 0\\0 & 2\end{array}\right] 
# \textrm{ et }
# {\bf (A^T A)^{-1}}=\left [\begin{array}{cc}1/3 & 0\\0 & 1/2\\\end{array}\right]
# $
# la solution aux moindres carrés est donc $\alpha=6$ et $\beta=5/2$.
# Le vecteur des erreurs ${\bf Y}-{\bf Ax}=(1/2$ $-1$ $1/2)^T $ est bien orthogonal dans $R^3$ aux colonnes de ${\bf A}$. Chaque erreur peut être également représentée par la distance verticale entre la mesure ${ Y_i}$ et la droite $f(t)=6+\frac{5}{2}t$ pour $t=t_i$ (segments rouges).
# ```{margin} 
# ![](./images/projection.png)
# ```
# ````
# 
# 
# 
# ````{prf:example}
#  Le système $\mathcal{S}$ produit pour les entrées $t_i =\{-2,-1,0,1,2\}$ les sorties \\$Y_i=\{3,1,-0.5,1.5,4\}$.
# \\On modélise $\mathcal{S}$ par un polynôme $P(t)=\alpha+\beta t+\gamma t^2 = {\bf a^T}{\bf x}$, où ${\bf x}=\begin{pmatrix}\alpha\\\beta\\\gamma\end{pmatrix} $et ${\bf a}=\begin{pmatrix}1\\t\\t^2\end{pmatrix}$. On recherche les valeurs de $\alpha,\beta,\gamma$ qui minimisent la somme des carrés des erreurs entre les valeurs théoriques et les valeurs mesurées.
# \begin{center}
# \begin{tabular}{|c|c|c|c|}
# \hline
# Point&$t_i$&$Y_i$&$y_i$\\
# \hline
# $P_1$&-2&3&$\alpha-2\beta+4\gamma$\\
# $P_2$&-1&1&$\alpha-\beta+\gamma$\\
# $P_3$&0&-0.5&$\alpha$\\
# $P_4$&1&1.5&$\alpha+\beta+\gamma$\\
# $P_5$&2&4&$\alpha+2\beta+4\gamma$\\
# \hline
# \end{tabular}
# \end{center}
# En posant ${\bf A}=\begin{pmatrix}
# 1&-2&4\\1&-1&1\\1&0&0\\1&1&1\\1&2&4
# \end{pmatrix}$ et ${\bf Y}\begin{pmatrix}
# 3\\1\\-0.5\\1.5\\4
# \end{pmatrix}
# $, on cherche ${\bf x}$ qui minimise $\norme{{\bf Ax}-{\bf Y}}^2$, soit encore ${\bf x}$ solution du système aux équations normales ${\bf A^T Ax}={\bf A^T Y}$\\
# On a alors $${\bf A^T A}=\begin{pmatrix}5&0&10\\0&10&0\\10&0&34\end{pmatrix}\quad{\bf A^T Y}=\begin{pmatrix}9\\2.5\\30.5\end{pmatrix}\textrm{ d'où } P(t)=\frac{1}{70}+\frac{1}{4}t+\frac{25}{28}t^2$$
# ````
# 
# ## Systèmes incompatibles
# ```{index} Système ; incompatible
# ```
# Soit un système linéaire incompatible ${\bf Ax}={\bf b}$, où ${\bf A}\in\mathcal{M}_{mn}(\mathbb R)$ est telle que $\rang({\bf A})<m$ et ${\bf b}\notin Im({\bf A})$. On supposera par exemple (comme dans le cas des moindres carrés) que $m>n$ et $\rang({\bf A})=n$. Le système n'a donc pas de solution, et on le remplace par le système aux équations normales obtenu en le multipliant par ${\bf A^T} $ : 
# $${\bf A^T Ax}={\bf A^T b}$$
# Ce système est en général mal conditionné car $\sigma({\bf A^T A}) = \sigma^2({\bf A})$ dans le cas d'une matrice carrée ${\bf A}$. La méthode de Gauss risque d'être inefficace et on lui préférera des méthodes basées sur des transformations orthogonales qui ont l'avantage d'être numériquement stables.\\\\
# Finalement, si les colonnes de ${\bf A}$ sont orthonormées ({\it i.e.} orthogonales deux à deux et de norme 1) ${\bf A^T A}=\mathbb I$ (l'identité dans la \og petite\fg{} dimension $n$) et la solution des équations normales est simplement ${\bf x^*}={\bf A^T b}$. La  \og bonne\fg{} stratégie pour résoudre le problème des moindres carrés est donc de construire une base orthonormée de $Im({\bf A})$ pour calculer explicitement la projection. On verra au paragraphe \ref{trortho} que cette construction revient à triangulariser la matrice par des transformations orthogonales.
