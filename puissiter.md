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

```{code-cell} ipython3
try:
    import sympy 
except ModuleNotFoundError: 
    !pip3 install --quiet sympy
    import sympy
import numpy, matplotlib
```

# Méthode des puissances itérées


## Quotient de Rayleigh

```{margin} 
![](./images/rayleigh.png)
```

```{index} Rayleigh;quotient de 
```
```{prf:definition}
Pour une matrice symétrique $A$, le *quotient de Rayleigh* est le rapport défini pour tout $x\neq 0$ par :

$\rho_A(x)=\frac{x^\top Ax}{x^\top x}$
```
On vérifie immédiatement que si $x$ est vecteur propre, le quotient de Rayleigh fournit la valeur propre associée : en effet $Ax=\lambda x \Rightarrow x^\top Ax=\lambda x^\top x$.

Si $\lambda_1$ et $\lambda_n$ sont respectivement la plus petite et la plus grande valeur propre de $A$, et $x^1$,$x^n$ les vecteurs propres associés, on a également les résultats suivants :

```{prf:theorem} Théorème min-max de Courant-Fischer
$\begin{align*}
\lambda_1&=\rho_A(x^1)=\displaystyle\min_{x\in \mathbb{R}^n}\{\rho_A(x)\}\\
\lambda_n&=\rho_A(x^n)=\displaystyle\max_{x\in \mathbb{R}^n}\{\rho_A(x)\}
\end{align*}
$

De plus, si les valeurs propres sont rangées dans l'ordre croissant, on a 

$\begin{align*}
\lambda_i&=\displaystyle\min_{S_i}\{\displaystyle\max_{x\in S_i}\{\rho_A(x)\}\}\\
\lambda_i&=\displaystyle\max_{S_{i-1}}\{\displaystyle\min_{x\in S_{i-1}^\bot}\{\rho_A(x)\}\}
\end{align*}$

où $S_i$ est un sous-espace quelconque de dimension $i$.
```

Le sous-espace $S_i$ pour lequel le quotient de Rayleigh est maximum est le sous-espace propre associé aux $i$ premières valeurs propres. Le sous-espace pour lequel il est minimum est orthogonal au sous-espace propre associé aux $i-1$ premières valeurs propres. C'est donc le sous-espace engendré par les $n-i+1$ vecteurs propres associés à $\{\lambda_i\cdots\lambda_n\}$.

## Méthode des puissances itérées

La méthode des puissances itérées permet de calculer le vecteur propre associé à la plus grande valeur propre.

Supposons $A$ symétrique de valeurs propres ordonnées selon $|\lambda_1|\leq\cdots|\lambda_{n-1}|<|\lambda_n|$.

On considère l'itération suivante définie à partir d'un vecteur initial $q_0$ donné, tel que $\|q_0\|=1$, et $q_0$ n'est pas orthogonal à $v^n$, le vecteur propre associé à la plus grande valeur propre isolée $\lambda_n$ :

$\begin{align*}
x_{k+1}&=Aq_k\\
q_{k+1}&=\frac{x_{k+1}}{\|x_{k+1}\|}
\end{align*}$

Par récurrence, on montre alors que  $q_k=\frac{A^kq_0}{\|A^kq_0\|}$ et comme les vecteurs propres $\{v^1\cdots v^n\}$ forment une base de $\mathbb{R}^n$, on peut écrire 

 $q_0=\displaystyle\sum{i=1}^n\alpha_iv^i,\quad\alpha_n\neq 0$
 et 

 $A^kq_0=\alpha_n\lambda_n^k\left (v^n+\displaystyle\sum_{i=1}^{n-1}\frac{\alpha_i}{\alpha_n}\left (\frac{\lambda_i}{\lambda_n}\right )^kv^i\right )$

 Lorsque $k\rightarrow\infty$, les rapports $\left (\frac{\lambda_i}{\lambda_n}\right )^k$ tendent vers 0 pour $i\neq n$, 
 ce qui signifie que la suite des itérés $\{q_k\}$ converge vers le vecteur propre $v^n$ ou $-v^n$. 
 On peut montrer de plus que $\|Aq_k\|$ tend vers $|\lambda_n|$ et que la convergence est linéaire de taux $\left |\frac{\lambda_{n-1}}{\lambda_n}\right |$ si $\alpha_{n-1}\neq 0$.




```{code-cell} ipython3
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from matplotlib.text import Annotation


class Arrow3D(FancyArrowPatch):

    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)
        
    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs) 

def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)
setattr(Axes3D, 'arrow3D', _arrow3D)


def plot_vector3d(ax,vector3d, origin=[0, 0, 0], **options):
    return ax.arrow3D(origin[0], origin[1], origin[2], vector3d[0], vector3d[1],vector3d[2],
              arrowstyle="-|>",
              **options)


class Annotation3D(Annotation):

    def __init__(self, text, xyz, *args, **kwargs):
        super().__init__(text, xy=(0, 0), *args, **kwargs)
        self._xyz = xyz

    def draw(self, renderer):
        x2, y2, z2 = proj_transform(*self._xyz, self.axes.M)
        self.xy = (x2, y2)
        super().draw(renderer)
        
        
def _annotate3D(ax, text, xyz, *args, **kwargs):
    '''Add anotation `text` to an `Axes3d` instance.'''

    annotation = Annotation3D(text, xyz, *args, **kwargs)
    ax.add_artist(annotation)
setattr(Axes3D, 'annotate3D', _annotate3D)

def puissiter(A,v0,lam,niter,epsilon):
    v = v0
    vv=[v0]
    l = np.dot(v0,np.dot(A,v0))
    ll = [l]
    k=0
    while np.fabs(lam-l)>epsilon and k<niter:
        w = np.dot(A,v)
        v = w/np.linalg.norm(w)
        l = np.dot(v,np.dot(A,v))
        vv.append(v)
        ll.append(l)
        k=k+1
    return ll, vv,k

A = np.array([[2.,1,-1],[1,3,1],[-1,1,4]])
w,v=np.linalg.eig(A)
vmax = v[:, np.argmax(w)]
lam =np.max(w)
print("La plus grande valeur propre de A est ",lam)

epsilon = 1e-4
niter=50
ll, vv,k = puissiter(A,np.ones(3),lam,niter,epsilon)

fig = plt.figure(figsize=(5,5))
plt.plot(range(len(ll)),ll,'-o',label='Puissances itérées')
plt.plot(range(len(ll)),lam*np.ones((len(ll)), dtype=np.uint8) ,'r')
plt.ylabel('valeur propre')
plt.xlabel('Iteration');
plt.text(0, lam+0.3, "$\lambda$", color="r", fontsize=18)
plt.legend()
plt.title("Valeur propre approchée à "+ str(epsilon)+" près en "+str(k)+" itérations")

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111, projection='3d')
#plt.subplot(122)
plot_vector3d(ax,vv[0], color="b",mutation_scale=20)
ax.annotate3D('$v_0$', (vv[0][0],vv[0][1],vv[0][2]), xytext=(10, 10), textcoords='offset points',color="b",fontsize=18)
for i in range(1,k):
    plot_vector3d(ax,vv[i], color="g",alpha = 1-float(i)/k,mutation_scale=20)
    if i%2==0 :
        ax.annotate3D('v('+str(i)+')', (vv[i][0],vv[i][1],vv[i][2]), xytext=(3, 3), textcoords='offset points')

#        plt.text(vv[i][0],vv[i][1]+0.05,vv[i][2],'v('+str(i)+')',color="g")
plot_vector3d(ax,vmax, color="r",mutation_scale=20)
ax.annotate3D('$v_\lambda$', (vmax[0],vmax[1],vmax[2]), xytext=(-8, 10), textcoords='offset points',color="r",fontsize=18)
plt.title("Vecteur propre approché")
plt.tight_layout()
```



```{prf:example} Suite de Fibonacci
Considérons la suite de  Fibonacci : 0, 1, 1, 2, 3, 5, 8, 13, ...
On cherche, sans calculer explicitement le terme, quelle va être la valeur du 100-ième terme 
La suite est définie par $ {F}_{k+2}={F}_{k+1}+{F}_{k} $, avec 
$F_0=0, F_1=1$.  En ajoutant l'équation  $F_{k+1}=F_{k+1}$, on note alors  $u_{k}$ 

${u}_{k}=\begin{bmatrix} {F}_{k+1} \\ {F}_{k} \end{bmatrix} $
Ainsi

$ {u}_{k+1}=\begin{bmatrix} 1 & 1 \\ 1 & 0 \end{bmatrix} \begin{bmatrix} {F}_{k+1} \\ {F}_{k} \end{bmatrix}=\begin{bmatrix} {F}_{k+1}+{F}_{k} \\ {F}_{k+1} \end{bmatrix} \\ {u}_{k+1}=\begin{bmatrix} 1 & 1 \\ 1 & 0 \end{bmatrix}{u}_{k} $
```

```{code-cell} ipython3
from sympy import Matrix

A = Matrix([[1, 1], [1, 0]])
A.eigenvals()
```


```{code-cell} ipython3
A.eigenvects()
```

```{code-cell} ipython3
S, D = A.diagonalize()
```

Et puisque ${u}_{k} = {A}^{k}{u}_{0}$ 
et $u_0=\begin{bmatrix}1\\0\end{bmatrix}$


```{code-cell} ipython3
u_zero =  Matrix([1, 0])
u_100 = A ** 100 * u_zero
u_100
```

## Méthode des puissances inverses

Pour les mêmes raisons, l'itération 

$\begin{align*}
Ax_{k+1}&=q_k\\
q_{k+1}&=\frac{x_{k+1}}{\|x_{k+1}\|}
\end{align*}$

avec $\|q_0\|=1$, et $q_0$ n'est pas orthogonal à $v^1$, converge vers la direction du vecteur propre associé à la plus 
petite valeur propre en module. On remarquera le coût de calcul en $O(n^2)$.


```{code-cell} ipython3
import matplotlib.pyplot as plt

def puissinverse(A,v0,lam,niter,epsilon):
    v = v0
    I = np.eye(len(v0))
    vv = [v0]
    l = np.dot(v0,np.dot(A,v0))
    ll = [l]
    k=0
    while np.fabs(lam-l)>epsilon and k<niter:
        w = np.linalg.solve(A-lam*I,v)
        v = w/np.linalg.norm(w)
        l = np.dot(v,np.dot(A,v))
        vv.append(v)
        ll.append(l)
        k=k+1
    return ll, vv,k
```

```{prf:remark}
:class: dropdown
1. Accélération par *décalage* : la matrice $A+\alpha I$ a les mêmes vecteurs propres que $A$ et ses valeurs propres sont décalées de la quantité $\alpha$. La méthode des puissances itérées inverses converge d'autant plus vite que les rapports $\left|\frac{\lambda_1}{\lambda_2}\right|^k$ tendent rapidement vers 0. On a donc intérêt à ce que $\lambda_1$ soit le plus proche possible de 0, et de plus, la méthode sera d'autant plus rapide que l'écart entre les deux plus petites valeurs propres se creuse. La technique du décalage consiste donc à remplacer $A$ par $A+\alpha I$, avec $\alpha\approx -\lambda_1$. Plus l'estimation de $\lambda_1$ sera précise, plus la convergence sera rapide. Toutefois, il faut que $\alpha\neq -\lambda_1$ pour éviter que la matrice ne devienne singulière
2. Technique de *déflation* : la méthode des puissances itérées peut être étendue pour permettre le calcul de toutes les valeurs propres d'une matrice symétrique. Supposons en effet calculée la plus grande valeur propre $\lambda_n$ ainsi qu'un vecteur propre associé $v^n$. Soit $P_n$ la la matrice de projection orthogonale sur l'hyperplan $(v^n)^\bot$. La matrice $P_nA$ possède les mêmes vecteurs propres que $A$ et les mêmes valeurs propres à l'exception de $\lambda_n$ qui est remplacée par 0. L'application de la méthode des puissances itérées à $P_nA$ permettra donc de calculer la deuxième plus grande valeur propre de $A$. Cette technique, dite de déflation, permet théoriquement de calculer toutes les valeurs propres de $A$. Elle est toutefois numériquement instable sans précautions, et on lui préférera généralement la méthode des puissances groupées.
```

