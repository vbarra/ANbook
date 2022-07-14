#!/usr/bin/env python
# coding: utf-8

# # Introduction
# 
# ## Introduction
# L'analyse spectrale est l'étude des valeurs propres et des vecteurs propres d'une matrice carrée. \\
# Les valeurs propres\index{valeur propre}
# d'une matrice carrée $A (n\times n)$ sont les $n$ solutions dans $\mathbb{C}$ de l'équation caractéristique
# 
# $det(\lambda I-A) = 0$
# 
# Du point de vue de l'algèbre linéaire, cela signifie que le noyau de $\lambda I-A$ contient des vecteurs non nuls, appelés vecteurs propres\index{vecteurs!propres} associés à $\lambda$.
# Donc, si $x\neq 0$ est un vecteur propre  de $A$ associé à la valeur propre $\lambda$,
# on a $Ax=\lambda x$.
# 
# Le calcul des $n$ solutions de l'équation caractéristique est coûteux dès que $n>2$ et le théorème d'Abel montre qu'on ne peut 
# espérer la résoudre par des radicaux
# dès que $n>4$. On recherchera donc des méthodes itératives qui permettent d'approcher ces racines et non de les calculer explicitement car, à la différence des méthodes de résolution de systèmes linéaires vues dans le chapitre 2, la convergence sera ici asymptotique. En fait, les méthodes qui seront présentées pour le calcul des valeurs propres sont utilisées pour extraire les racines d'un polynôme en passant par la matrice compagne :  
# $$\displaystyle\sum_{i=0}^{n-1}a_it^i+t^n = 0$$
# est le polynôme caractéristique
# de la matrice compagne
# 
# $A =
# \left (
# \begin{array}{llll}
# 0 & \cdots &  &-a_0\\
# 1  &\ddots && -a_1\\
# 0 & & 0 &\vdots\\
# 0 & & 1 & -a_{n-1}
# \end{array}
# \right )
# $
# ```{index} Matrice;compagne
# ```
# ```{index} Polynôme;caractéristique
# ```
# 
# ## Intérêts
# 
# Considérons par exemple un système d'équations différentielles ordinaires : \index{système!différentiel}\\
# Trouver les fonctions $v(t)$ et $w(t)$, pour $t\in \mathbb{R}^+$, telles que :
# \begin{eqnarray*}
# \frac{dv}{dt}& =& 4v - 5w\\[.25pc]
# \frac{dw}{dt}& =& 2v - 3w
# \end{eqnarray*}
# avec comme conditions initiales $v(0)=8,w(0)=5$.
# Ce système peut être mis sous la forme vectorielle
# \begin{equation}\label{eqvec}
# \frac{du}{dt}=Au
# \end{equation}
# avec : 
# $u(t)=\begin{pmatrix}v(t)\\w(t)\end {pmatrix},u(0)=\begin{pmatrix}8\\5\end {pmatrix},A=\begin{pmatrix}4\quad -5\\2\quad -3\end {pmatrix}$
# 
# Sachant que la solution de l'équation $x'=ax$ en dimension 1 est $x(t)=x_0e^{at}$, qui diverge si $a>0$ et se stabilise asymptotiquement à zéro si $a<0$, on cherchera des solutions particulières de la forme 
# 
# $
# \left \{
# \begin{array}{cc}
# v(t) =& \tilde{v}e^{\lambda t}\\
# w(t) =& \tilde{w}e^{\lambda t}
# \end{array}
# \right .
# $ 
# En remplaçant dans (\ref{eqvec}), on voit que 
# 
# $
# u(t)=e^{\lambda t}\begin{pmatrix}
#                             \tilde{v}\\ \tilde{w}
#                   \end{pmatrix}
#     \Rightarrow \frac{du}{dt}=\lambda e^{\lambda t}
#                   \begin{pmatrix}
#                              \tilde{v}\\ \tilde{w}
#                   \end{pmatrix}=Au(t)
# $
# d'où 
# $
# A\begin{pmatrix}\tilde{v}\\\tilde{w}\end{pmatrix}
# =\lambda\begin{pmatrix}\tilde{v}\\\tilde{w}\end{pmatrix}
# $
# 
# et $\lambda$ doit donc être une valeur propre de $A$ associée au vecteur propre $\begin{pmatrix}\tilde{v}\\\tilde{w}\end {pmatrix}$. Dans le cas présent, on calcule aisément les éléments propres de $A$ :
# 
# $\begin{align*}
# \lambda_1=-1 & \Rightarrow \begin{pmatrix}\tilde{v_1}\\\tilde{w_1}\end {pmatrix}=\begin{pmatrix}1\\1\end {pmatrix} \Rightarrow u^1(t)=\begin{pmatrix}e^{-t}\\e^{-t}\end {pmatrix}\\
# \lambda_1=2 & \Rightarrow \begin{pmatrix}\tilde{v_2}\\\tilde{w_2}\end {pmatrix}=\begin{pmatrix}5\\2\end {pmatrix} \Rightarrow u^2(t)=\begin{pmatrix}5e^{2t}\\2e^{2t}\end {pmatrix}
# \end{align*}$
# 
# La linéarité du système implique que toute combinaison linéaire des solutions particulières $u^1,u^2$ est solution de l'équation différentielle. La solution générale s'écrit donc :
# 
# $u(t)=\alpha u^1(t)+\beta u^2(t)$
# 
# où $\alpha=3$ et $\beta=1$ sont déterminées par les conditions initiales.
