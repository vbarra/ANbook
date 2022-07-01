#!/usr/bin/env python
# coding: utf-8

# # Tutoriel Matplotlib

# [Matplotlib](https://matplotlib.org) est une librairie d'affichage graphique. On ne s'intéresse ici qu'au module `matplotlib.pyplot`

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# on peut afficher une courbe à partir d'une liste de points, d'une ou plusieurs fonctions, faire varier le type et la couleur des courbes, ajouter des légendes, des titres, la sauver dans un fichier

# In[2]:


plt.plot([1, 2, 4, 9, 5, 3])
plt.title("Premier graphique")
plt.show()

x = np.arange(0, 3 * np.pi, 0.2)
y_sin = np.sin(x)
y_cos = np.cos(x)
y_sqr = np.square(x)
plt.plot(x, y_sin,"g--")
plt.plot(x, y_cos,'r-')
plt.plot(x,x**2/50, 'b^')

plt.xlabel('x axis label')
plt.ylabel('y axis label')
plt.legend(['Sinus', 'Cosinus','Carré'])

plt.title("Fonctions cosinus et sinus")
plt.savefig("plots.png", dpi=50)


# On peut utiliser une même figure pour afficher plusieurs graphiques côte à côte ([documentation](http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.subplot))

# In[3]:


x = np.linspace(-1.4, 1.4, 30)
plt.subplot(2, 2, 1)  
plt.title('y=x')
plt.plot(x, x)
plt.subplot(2, 2, 2)  
plt.plot(x, x**2)
plt.title('y=x*x')
plt.subplot(2, 2, 3)  
plt.title('y=x^3')
plt.plot(x, x**3)
plt.subplot(2, 2, 4)  
plt.title('y=x^4')
plt.plot(x, x**4)
plt.tight_layout()


# In[ ]:




