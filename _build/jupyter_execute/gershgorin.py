#!/usr/bin/env python
# coding: utf-8

# # Disques de Gershgorin
# 
# Avant d'aborder quelques algorithmes de calcul des valeurs propres d'une matrice, donnons une alternative pratique à ces algorithmes . Le théorème suivant  permet de localiser les valeurs propres dans des disques, dits disques de Gershgorin, du plan complexe.
# 
# ```{prf:theorem} Théorème de Gershgorin
# Si on représente une matrice $A$ (ou toute matrice semblable à $A$) sous la forme $A=diag\{d_1\cdots d_n\}+F$, où $F$ est une matrice de diagonale nulle, alors le spectre de $A$ est contenu dans l'union des disques $D_i,1\leq i\leq n$ du plan complexe, tels que 
# 
# $D_i=\left \{ z\in \mathbb{C}, |z-d_i|\leq \displaystyle\sum_{j=1}^n|f_{ij}|\right \}$
# ```
# 
# ```{index} Gershgorin;théorème
# ```
# 
# Une application intéressante de ce résultat est l'estimation des valeurs propres d'une matrice obtenue en perturbant une matrice dont on connaît le spectre.
# 
# 
# ```{margin} 
# ![](./images/gershgorin.png)
# ```
# 
# ```{prf:example}
# $A =
# \left[
# \begin{array}{ccc}
# 1&0.1&-0.1\\
# 0&2&0.4\\
# -0.2&0&3\\
# \end{array}
# \right]
# $
# 
# dont les valeurs propres sont situées dans les disques suivants
# $\begin{align*}
# D_1&=\left \{ z\in \mathbb{C}, |z-1|\leq 0.2\right\}\\
# D_2&=\left \{ z\in \mathbb{C}, |z-2|\leq 0.4\right\}\\
# D_3&=\left \{ z\in \mathbb{C}, |z-3|\leq 0.2\right\}
# \end{align*}$
# ```

# In[1]:


import numpy as np
N = 5 # dimension of our square matrix

D = np.diag([0, 3 + 1j, 4 + 1j, 1 + 5j, 9 + 2j])
M = np.random.rand(N, N) + D

R = np.zeros(N) # disk radii
for i in range(N):
    R[i] = sum(abs(M[i,:])) - abs(M[i,i])

eigenvalues = np.linalg.eigvals(M)

# Plotting code
fig, ax = plt.subplots()
for k in range(N):
    x, y = M[k,k].real, M[k,k].imag
    ax.add_artist( plt.Circle((x, y), R[k], alpha=0.5) )
    plt.plot(eigenvalues[k].real, eigenvalues[k].imag, 'k+')

ax.axis([-4, 12.5, -4, 9])
ax.set_aspect(1)    
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("Gershgorin disks and eigenvalues $x + iy$")

