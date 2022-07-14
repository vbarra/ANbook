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

# Résultats généraux

On donnera ici quelques propriétés d'intérêt surtout pratique concernant les valeurs propres et les vecteurs propres de certaines matrices. Les démonstrations sont omises.

Une matrice carrée $A$ de dimension $n$ à coefficients complexes ou réels possède $n$ valeurs 
propres non nécessairement distinctes dans $\mathbb{C}$. L'ensemble de ces valeurs propres est 
le *spectre* de $A$, noté $\mathrm{sp}(A)$.
```{index} Matrice;spectre
```
```{index} Spectre
```
 Le nombre de fois où apparaît une 
valeur propre $\lambda$ dans ce spectre est appelé sa 
\textbf{multiplicité}.
```{index} Valeur propre;multiplicité
```
```{index} Rayon spectral
```
Le *rayon spectral*, noté $\rho(A)$, est le plus grand module des valeurs propres de $A$. La somme des valeurs propres est égale à la \textbf{trace} de la matrice :

$\displaystyle\sum_{i=1}^n\lambda_i=\displaystyle\sum{i=1}^na_{ii}$

et on en déduit donc que si une matrice réelle possède une valeur propre complexe, son conjugué est aussi valeur propre. De même, le produit des valeurs propres est égal au \textbf{déterminant} de 
$A$ :
$\displaystyle\prod_{i=1}^n\lambda_i=det(A)$

Donc $A$ non singulière $\Leftrightarrow$ $0\notin \mathrm{sp}(A)$.

```{warning}
Il n'est pas possible en général d'obtenir valeurs et vecteurs propres de la somme ou du produit de deux matrices dont on connaît le spectre, mais certains cas particuliers où les vecteurs propres sont les mêmes permettent de conclure :

- si $\lambda$ est valeur propre de $A$,  $\lambda^k$ est valeur propre de $A^k$ (en particulier, si $A$ est inversible, $\lambda^{-1}$ est valeur propre de l'inverse)
- si $\lambda$ est valeur propre de $A$, $\lambda+\alpha$ est valeur propre de $A+\alpha I$.

Si $\lambda$ est une valeur propre de $A$, le système $(A-\lambda I)x=0$ possède des solutions non nulles appelées \textbf{vecteurs propres} 
de $A$ associés à $\lambda$. On appelle \textbf{sous-espace propre}
associé à une valeur propre $\lambda$ le noyau de $A-\lambda I$. Sa dimension est donc au moins égale à 1.
```
```{index} Vecteur propre
```
```{index} Sous-espace;propre
```

Lorsque les valeurs propres sont distinctes, les vecteurs propres sont linéairement indépendants. L'implication réciproque est fausse, on pensera au cas de la matrice identité pour s'en 
convaincre.

Si les $n$ vecteurs propres $x_i,1\leq i\leq n$ peuvent être choisis linéairement indépendants, ils forment une matrice $X=\left [ x_i\right ]$ non singulière, et on a :

$X^{-1}AX=\Lambda=\mathrm{diag}\{\lambda_1\cdots\lambda_n\}$

En effet, $AX$ est la matrice dont les colonnes sont les vecteurs $Ax_i=\lambda x_i$ qui est bien égale à $X\Lambda$. On dit dans ce cas que $A$ est *diagonalisable*
```{index} Matrice;diagonalisable
```
Cette propriété ne peut s'écrire qu'avec la matrice des vecteurs propres. 
Elle caractérise le fait que $X$ a pour colonnes des vecteurs propres et que ces vecteurs 
propres sont linéairement indépendants.


```{index} Matrice;défective
```
```{admonition} Observation
Il existe des matrices qui ne sont pas diagonalisables (on les appelle 
matrices *défectives*. 

Elles satisfont aux deux conditions 
suivantes :

- il existe des valeurs propres multiples
- la dimension du sous-espace propre associé est strictement inférieure à la multiplicité de la valeur propre.
```

Le seul fait de l'existence de valeurs propres multiples ne suffit pas à impliquer que la matrice soit défective (cf. la matrice identité par exemple).
Un exemple typique de matrice défective est $A=\begin{pmatrix} 0\quad 1\\0\quad 0\end{pmatrix}$, qui possède la valeur propre double 0. Les vecteurs propres sont de la forme $x=\begin{pmatrix}\alpha\\0\end{pmatrix}$, ils engendrent donc un sous-espace de dimension 1. $A$ n'est donc pas diagonalisable.

Dans certains cas, il est relativement facile d'inférer sur la valeur ou la nature des valeurs propres d'une matrice $A$:
- si $A$ est diagonale, les valeurs propres sont les éléments diagonaux
- si $A$ est triangulaire, les valeurs propres sont également les éléments diagonaux
- si $A$ est non singulière, les valeurs propres sont toutes différentes de 0
- si $A$ est orthogonale, les valeurs propres ont pour module 1 et il est possible de choisir une base de vecteurs propres orthonormés
- si $A$ est symétrique, ses valeurs propres sont réelles. Les vecteurs propres associés à des valeurs propres distinctes sont 
alors orthogonaux. On a $\|A\|_2=\rho(A)$.

Le dernier point implique en particulier que toute matrice réelle symétrique est diagonalisable par une matrice orthogonale : $\Lambda=Q^\top AQ$,où $Q$ est formée par $n$ vecteurs propres orthonormés. On a alors la *factorisation spectrale* 
d'une matrice symétrique réelle :
$A=Q\Lambda Q^\top $
```{index} Factorisation spectrale
```
```{important}
En résumé de ce qui précède :
- l'inversibilité est liée aux valeurs propres
- la diagonalisabilité est liée aux vecteurs propres
```
