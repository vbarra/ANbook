#!/usr/bin/env python
# coding: utf-8

# # Tutoriel Numpy

# [Numpy](https://numpy.org) est la bibliothèque centrale de calcul scientifique en Python. Elle utilise un type de tableau multidimensionnel, et des outils pour travailler sur ces tableaux.
# Pour utiliser la bibliothèque `numpy`, on l'importe :

# In[1]:


import numpy as np


# ## Arrays
# ### Définition
# Un `array` est un tableau de valeurs, toutes du même type, indexées par un tuples d'entiers. Le nombre de dimensions est le rang du tableau. 

# In[2]:


# Création d'un vecteur
v = np.array([1, 2, 3])  
print('taille du vecteur : ',v.shape)
for i in range(v.shape[0]):
    print("l'élément ",i," du vecteur est :",v[i])

# Matrice      
m = np.array([[1,2,3],[4,5,6]])   
print('\ntaille de la matrice : ',m.shape)
print(m)
     


# In[3]:


# Initialisations particulières
z = np.zeros((2,2))         # matrice nulle
o = np.ones((5,3))          # matrice de 1
c = np.full((2,2), 7)       # matrice de 7
i = np.eye(4,4)              # matrice identité
r = np.random.random((2,2)) # valeurs aléatoires


# ### Indexation
# On résume dans la suite quelques méthodes d'accès, voir la [documentation](https://docs.scipy.org/doc/numpy/reference/arrays.indexing.html) pour plus de détails.

# * par Slicing : 

# In[4]:


a = np.array([[1,2,3,4], [5,6,7,8], [9,10,11,12]])
print(a)
b = a[:2, 1:3]
print("Sous-tableau composé des coefficients à l'intersection  deux premières lignes et des colonnes 1 et 2\n", b)
print('a[0, 1] = ',a[0, 1])  
b[0, 0] = -100    
print('b pointe sur a, donc modifier b modifie a :nouvelle valeur de a[0, 1] : ',a[0, 1]) 


# * par accès aux lignes/colonnes : 

# In[5]:


lig_r1 = a[1, :]    # deuxième ligne de a, de range 1
lig_r2 = a[1:2, :]  # deuxième ligne de a, de range 2
lig_r3 = a[[1], :]  # deuxième ligne de a, de range 2
print('ligne ',lig_r1,'de rang ', lig_r1.shape) 
print('ligne ',lig_r2, 'de rang ',lig_r2.shape)
print('ligne ',lig_r3, 'de rang ',lig_r2.shape)

col_r1 = a[:, 1]
col_r2 = a[:, 1:2]
print('colonne ',col_r1, 'de rang ',col_r1.shape)
print('colonne ',col_r2, 'de rang ',col_r2.shape)


# * par tableau d'indices

# In[6]:


a = np.array([[1,2,3], [4,5,6], [7,8,9], [10, 11, 12]])
print("a=",a)
# On sélectionne le 1e élément de la première ligne, le 3e de la deuxième ligne, le 1e de la troisième et le 2e de la quatrième
b = np.array([0, 2, 0, 1])
res =  a[np.arange(4), b]
print("Sélection par une liste d'indices", res)  

a[np.arange(4), b] -=100
print("On ne modifie que les coefficients décrits par b")
print(a)


# * par accès booléen

# In[7]:


a = np.array([[1,2,3], [4,5,6], [-1,-1,-3], [1, 2, 3]])
bool_idx = (a > 2)  
print(bool_idx)
# Accès aux éléments de a satisfaisant la contrainte
b = a[bool_idx]
c = a[a > 2]
print('résultat : ',b,' également accessible par ',c)


# ## Types de données
# `numpy` essaye de déterminer le type de données lors de la création, mais on peut spécifier explicitement le type désiré ([documentation](http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html)).

# In[8]:


x = np.array([1, 2])  
y = np.array([1.0, 2.0])  
z = np.array([1, 2], dtype=np.int64)  

print('type de x : ',x.dtype, '\ntype de y : ',y.dtype, '\ntype de z : ',z.dtype)


# ### Opérations mathématiques
# opèrent point à point dans le tableau. Sont disponibles comme surchare d'opérateur ou comme fonction ([documentation](http://docs.scipy.org/doc/numpy/reference/routines.math.html)).

# In[9]:


x = np.array([[1,2],[3,4]], dtype=np.float64)
y = np.array([[5,6],[7,8]], dtype=np.float64)

print('\n Addition')
print(x + y)
print(np.add(x, y))

print('\n Division')
print(x / y)
print(np.divide(x, y))

print('\n Racine carrée')
print(np.sqrt(x))


# On récupère également les opérations d'algèbre linéaire

# In[10]:


x = np.array([[1,2],[3,4]])
y = np.array([[5,6],[7,8]])

v = np.array([9,10])
w = np.array([11, 12])

# Transposée
print('\n Transposée')
print(x.T)

# Produit scalaire
print('\n Produit scalaire')
print(v.dot(w), ' ou ', np.dot(v, w))

#Produit matrice/vecteur
print('\n Produit matrice/vecteur')
print(x.dot(v), ' ou ', np.dot(x, v))

# Produit matrice/matrice
print('\n Produit matrice/matrice')
print(x.dot(y), ' ou \n', np.dot(x, y))


# D'autres fonctions sont classiquement utilisées comme par exemple `sum`:

# In[11]:


x = np.array([[1,2],[3,4]])

print("Somme de tous les éléments d'un tableau : ",np.sum(x))  
print("Somme de tous les éléments de chaque colonne d'un tableau : ",np.sum(x, axis=0)) 
print("Somme de tous les éléments de chaque ligne d'un tableau : ",np.sum(x, axis=1))  


#     
# ### Broadcasting 
# [Le broadcasting](http://docs.scipy.org/doc/numpy/user/basics.broadcasting.html) : permet de travailler avec des tableaux de différentes tailles pour réaliser des opérations arithmétiques. 
# 

# Par exemple, pour ajouter un vecteur constant à chaque ligne d'une matrice 

# In[12]:


x = np.array([[1,2,3], [4,5,6], [7,8,9], [10, 11, 12]])
v = np.array([1, 0, 1])
y = x + v  
print(y)

