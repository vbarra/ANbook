#!/usr/bin/env python
# coding: utf-8

# # Puissances groupées et méthode QR
# 
# On peut généraliser la méthode des puissances itérées à des itérations dites groupées, permettant d'identifier les vecteurs propres associés aux $p$ plus grandes valeurs propres
# en valeur absolue. A chaque itération, on applique la transformation $A$ à chacun des vecteurs orthonormés, et on orthogonalise (par Gram-Schmidt ou Householder) le système résultant.
# 
# ## Itération des puissances itérées
# 
# Soit $\{q_1^{(k)}\cdots q_p^{(k)}\}$ un système de $p$ vecteurs orthonormés obtenu à l'itération $k$ On calcule les $p$ premiers vecteurs propres en itérant les étapes suivantes : 
# 
# 1. construire le système $\{Aq_1^{(k)}\cdots Aq_p^{(k)}\}$
# 2. calculer les quotients de Rayleigh : $\lambda_i={q^{(k)}_i}^\top Aq_i^{(k)},\quad 1\leq i\leq p$
#     \item tester la convergence : $\displaystyle\max_{1\leq i\leq p} \|Aq_i^{(k)}-\lambda_i q_i^{(k)} \|<Tol$
# 3. orthogonaliser $\{Aq_1^{(k)}\cdots Aq_p^{(k)}\} \rightarrow \{q_1^{(k+1)}\cdots q_p^{(k+1)}\}$
# 
# 
# 
# ## Méthode QR
# 
# 
# La méthode QR correspond au cas $p=n$ : soit $Q_k$ la matrice ($n\times n$) dont les colonnes sont les vecteurs orthogonaux $q_i^{(k)}$ et soit $A_k=Q_k^\top AQ_k$ la représentation de $A$ sur la base des $\{q_i^{(k)}\}$. On peut écrire alors (Gram-Schmidt) : 
# 
# $AQ_k=Q_{k+1}R_{k+1}$ 
# 
# où $R_{k+1}$ est triangulaire supérieure. Et ainsi 
# 
# $\begin{align*}
# A_k&=Q_k^\top Q_{k+1}R_{k+1}=QR_{k+1}\\
# A_{k+1}&=Q_{k+1}^\top AQ_{k+1}=R_{k+1}Q_k^\top Q_{k+1}=R_{k+1}Q
# \end{align*}$
# 
# En résumé, $Q$ est la matrice orthogonale qui permet de triangulariser $A_k$ et on écrit $A_{k+1}$ en inversant le produit $QR$. La méthode QR prend alors la forme particulièrement simple suivante :
# 
# ```{prf:algorithm} Algorithme QR
# :label: FactLU
# **Entrée : ** La matrice ${\bf A}$
# 
# **Sortie : ** Les matrices ${\bf Q} et {\bf R}
# 
# 1. $A_0=A$
# 2. Tant que le plus grand élément non diagonal de  $A_k$ est non nul
#     1. $A_k=Q_kR_k$ (par Gram-Schmidt par exemple)
#     2. $A_{k+1}=R_kQ_k$
#     
# ```
# 
# Observons que $A_{k+1}=Q^\top A_kQ$, ce qui implique que les matrices sont toutes semblables et que la suite des matrices  $Q$ converge vers la matrice des vecteurs propres. La matrice $A_k$ converge vers une matrice diagonale où les valeurs propres sont rangées dans l'ordre décroissant (Figure \ref{52}). 
# 
# 
# 
# ```{prf:remark}
# Les performances sont considérablement améliorées si on intègre deux modifications :
# - {\gr décalage} de la matrice en prenant comme approximation de la plus petite valeur propre l'élément ($n\times n$) de $A_k$
# - transformation préalable de $A$ en une matrice {\gr tridiagonale}. Les matrices $A_k$ restent tridiagonales et l'orthogonalisation s'effectue en $O(n)$ flops.
# ```
