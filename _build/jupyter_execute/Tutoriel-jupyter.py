#!/usr/bin/env python
# coding: utf-8

# # Tutoriel Jupyter

#  [Jupyter](https://jupyter.org/) est une application web qui permet aux utilisateurs de créer et de partager des documents appelés notebooks. Les notebooks Jupyter peuvent contenir du texte, des équations et des graphiques, ainsi que du code Python en direct avec lequel l'utilisateur peut interagir.

# ## Cellules
# Chaque notebook est constitué d'un certain nombre de cellules. Chaque cellule est soit une cellule de texte (en [Markdown](https://www.markdownguide.org)), soit une cellule de code. Les cellules markdown contiennent principalement du texte et des instructions [LaTex]. Lorsqu'elles sont exécutées, ces cellules génèrent l'ensemble du texte formaté et des formules. Les cellules de code contiennent des instructions, écrites (ici) en Python, qui doivent être exécutées par le Notebook.
# Toutes les cellules des notebooks, qu'elles soient de type markdown ou de type code, doivent être exécutées afin de produire des résultats. Pour exécuter une cellule, sélectionnez-la et appuyez sur "Shift + Enter". Vous pouvez sélectionner une cellule en cliquant dessus, ou en vous déplaçant vers le haut ou vers le bas avec les touches fléchées. Pour modifier une cellule, vous devez passer en mode édition. Lorsqu'il n'est pas en mode édition, on dit que le notebook est en mode commande.
# 
# Bien que toutes les cellules puissent être modifiées et exécutées, il existe certaines distinctions entre les deux types de cellules.
# - Une erreur Python s'affiche si le code Python n'est pas correct. Dans une cellule Markdown, il n'y a pas de message d'erreur, bien que le texte ou la formule affichés puissent sembler incorrects. Les cellules Markdown ne sont pas exécutées dans le même sens que les cellules de code.
# - Toutes les cellules Markdown sont exécutées lorsque le notebook est ouvert, mais pas les cellules de code, qui doivent être explicitement exécutées en appuyant sur 'Shift + Enter'. 
# 

# ## Menus
# Plusieurs menus déroulants sont disponibles dans un notebook.
# - Le menu Fichier contient des options familières permettant de créer ou d'ouvrir d'autres notebooks, d'enregistrer ou de renommer le notebook actuel, de télécharger ou de fermer le le notebook courant.
# - Le menu Editer contient plusieurs options permettant de manipuler les cellules du notebook. Les opérations courantes sur les cellules sont : couper, copier, ajouter, supprimer, fusionner, etc.
# - Le menu View contient quelques options d'affichage.
# - Le menu Cellule offre différentes options pour exécuter des collections de cellules, plutôt que des cellules individuelles. Il contient également l'option Cell Type, qui permet de basculer les cellules entre le code et le Markdown.
# - Le menu Kernel interagit avec l'interpréteur Python. Il est parfois nécessaire de redémarrer le noyau pour recharger un module, ou d'interrompre le noyau si qun problème survient, mais la plupart du temps, l'interaction avec le noyau n'est pas nécessaire pour les opérations de routine.
# - Le menu Widget interagit avec les widgets, qui sont des éléments dynamiques qui peuvent être ajoutés aux Notebooks. 
# - Le menu Aide offre de la documentation sur l'environnement du Notebook ainsi que des liens vers des références externes.
# 

# ## Jupyter book
# Le cours est proposé en format [Jupyterbook](https://jupyterbook.org/en/stable/intro.html), ce qui vous permet, en plus des notions théoriques, de voir des exemples de codes en Python et de les modifier (en cliquant sur la petite fusée en haut, les cellules de codes deviennent exécutables)

# ## Resources
# [A Gallery of Interesting Jupyter Notebooks](https://github.com/jupyter/jupyter/wiki/A-gallery-of-interesting-Jupyter-Notebooks)
# 
# [The Programming Historian](https://programminghistorian.org/en/lessons/jupyter-notebooks)
# 
# [Jupyter Notebook Documentation](https://jupyter-notebook.readthedocs.io/en/stable/notebook.html)
