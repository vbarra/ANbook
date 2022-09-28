#!/usr/bin/env python
# coding: utf-8

# # Transformations élémentaires
# Parmi toutes les méthodes de résolution de systèmes linéaires, nous étudions dans la suite deux algorithmes largement utilisés : l'élimination de Gauss et la méthode LU.
# 
# ## Définition
# ````{prf:definition} Transformation élémentaire
# ```{index} Transformation ; élémentaire
# ```
# Soit ${\bf A}\in\mathcal{M}_{mn}(\mathbb R)$. On appelle transformation élémentaire
# des lignes (respectivement des colonnes) de ${\bf A}$ une combinaison des trois transformations
# suivantes :
# - multiplication d'une ligne (resp. colonne) par un scalaire non nul;
# - remplacement d'une ligne(resp. colonne)  par la somme de cette ligne et d'une autre ligne
#   (resp. colonne);
# - permutation de deux lignes (resp. colonnes).
# ````
# 
# Une transformation élémentaire sur les lignes (resp. colonnes) de ${\bf A}$
# équivaut à multiplier ${\bf A}$ à gauche (resp. à droite) par une matrice
# élémentaire obtenue en appliquant la même transformation à la matrice
# identité (de taille $m$ pour les transformations sur les lignes, et de taille $n$ pour les transformations sur les colonnes). C'est une transformation non singulière, elle ne change donc pas le rang de ${\bf A}$.
# 
# 
# Si on note pour $\alpha\in\mathbb R^*$:
# 
# - ${\bf M_i(\alpha)}$ la matrice élémentaire qui réalise ${\bf A_{i,\bullet}}\leftarrow \alpha {\bf A_{i,\bullet}}$ (ou ${\bf A_{\bullet,i}}\leftarrow \alpha {\bf A_{\bullet,i}}$)
# - ${\bf C_{ij}(\alpha)}$ la matrice élémentaire qui réalise ${\bf A_{i,\bullet}}\leftarrow {\bf A_{i,\bullet}} + \alpha {\bf A_{j,\bullet}}$ (ou ${\bf A_{\bullet,i}}\leftarrow {\bf A_{\bullet,i}} + \alpha {\bf A_{\bullet,j}})$
# - ${\bf P_{ij}}$ la matrice élémentaire qui réalise ${\bf A_{i,\bullet}}\leftrightarrow {\bf A_{j,\bullet}}$  (ou ${\bf A_{\bullet,i}}\leftrightarrow {\bf A_{\bullet,j}}$)
# 
# Alors on a le résultat : 
# 
# 
# ```` {prf:theorem} 
# Les matrices ${\bf M_i(\alpha)}$, ${\bf C_{ij}(\alpha)}$ et ${\bf P_{ij}}$ sont inversibles et :
# 
# $({\bf M_i(\alpha))^{-1}} = {\bf M_i(\alpha^{-1})},\quad ({\bf C_{ij}(\alpha))^{-1}} =  {\bf C_{ij}(-\alpha)} \quad \textrm{et}\;   {\bf P_{ij}^{-1}}=  {\bf P_{ij}}$
# ````
# 
# ````{prf:example}
# Soit $
# {\bf A}=
# \begin{pmatrix}1 & 0 & 0 & 2\\ 1 & 2 & 0 & 1\\ 2 & 1 & 3 & 4\end{pmatrix}
# $
# 
# Si on additionne les deux premières lignes et que l'on inscrive le résultat
# sur la première ligne sans changer les deux autres, on obtient
# $
# {\bf A'}=
# \begin{pmatrix}2 & 2 & 0 & 3\\ 1 & 2 & 0 & 1\\ 2 & 1 & 3 & 4\end{pmatrix}
# $
# 
# On peut vérifier que ${\bf A'}={\bf C_{12}(1)}{\bf A}$, avec
# $
# {\bf C_{12}(1)}=\begin{pmatrix}1 & 1 & 0 \\ 0 & 1 & 0\\ 0 & 0 & 1\end{pmatrix}
# $
# 
# De même la permutation des première et dernière colonnes de ${\bf A}$ s'effectue en multipliant à droite ${\bf A}$ par 
# ${\bf P_{14}}=
# \begin{pmatrix}0&0&0 & 1  \\ 0 & 1 & 0&0\\ 0 & 0 & 1&0\\1&0&0&0\end{pmatrix}
# $
# ````
# 
# ## Pivotage d'une colonne
# ```{index} Pivotage
# ```
# 
# C'est une séquence de transformations élémentaires qui transforme une colonne
# d'une matrice ${\bf A}$ en vecteur unitaire. L'élément qui doit se transformer en 1 est le `pivot`.
# ```{index} Pivot
# ```
# L'algorithme {prf:ref}`pivotage` décrit la transformation d'une colonne ${\bf A_{\bullet s}}$ en le $s^e$ vecteur de la base canonique.
# 
# ```{prf:algorithm} Pivotage d'une colonne ${\bf A_{\bullet s}}$ suivant le pivot $a_{rs}\ne 0$
# :label: pivotage
# 
# **Entrée :** ${\bf A}$, $s$ l'indice de la colonne à pivoter, $r$ l'indice de ligne du pivot
# 
# **Sortie :** La colonne transformée
# 
# 1. ${\bf A_{r \bullet}}\leftarrow \frac{1}{a_{rs}}\bf A_{r \bullet}$
# 2. Pour $i \ne r$ : ${\bf A_{i \bullet}}\leftarrow {\bf A_{i \bullet}} -a_{is}{\bf A_{r \bullet}}$
# ```
# 
# 
# 
# Pour $m=n$, il est clair que, si tous les pivots successifs sont non nuls, une séquence de
# $n$ pivotages effectués sur les $n$ colonnes de la matrice avec le pivot
# sur la diagonale transformera la matrice en l'identité. Si ${\bf E_i}$ est la 
# matrice élémentaire associée au pivotage de la $i^e$ colonne, on a
# $
# \mathbb I_n={\bf E_n}{\bf E_{n-1}}\cdots {\bf E_1}{\bf A}
# $
# d'où une première méthode pour calculer l'inverse d'une matrice : 
# $
# {\bf A^{-1}}={\bf E_nE_{n-1}}\cdots {\bf E_1}.
# $
# 
# Donc pour calculer l'inverse d'une matrice, il suffit d'effectuer les pivotages
# en parallèle sur la matrice identité. La matrice à inverser se transforme
# progressivement en la matrice identité alors que l'identité devient 
# l'inverse. 
# 
# S'il est vrai que, quand aucun pivot nul n'est rencontré, cet
# algorithme a une complexité de $O(n^3)$, il peut atteindre $O(2n^3)$ dans le
# cas général et nécessite le stockage de deux matrices $n\times n$. On lui
# préférera les méthodes étudiées dans la suite, basées elles aussi sur
# des pivotages successifs de la matrice, mais plus robustes et moins 
# coûteuses.
# 
# ````{prf:example}
# Pivotons la première colonne de 
# $
# {\bf A}=
# \begin{pmatrix}1 & 1 & 1\\ 1 & 0 & -1\\ -2 & 1 & 3\end{pmatrix}
# $
# avec $a_{11}=1$ comme pivot. On doit donc multiplier $A$ à gauche par
# la matrice élémentaire 
# $
# {\bf E}=
# \begin{pmatrix}1 & 0 & 0\\ -1 & 1 & 0\\ 2 & 0 & 1\end{pmatrix}$
# et le résultat est 
# ${\bf A'}={\bf E}{\bf A}=
# \begin{pmatrix}1 & 1 & 1\\ 0 & -1 & -2\\ 0 & 3 & 5\end{pmatrix}
# $
# ````
