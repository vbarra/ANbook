import numpy as np 
import matplotlib.pyplot as plt
from numpy import linalg as LA

fig,ax= plt.subplots(3,2,figsize=(10,8))
nbpoints = 1000
c= np.linspace(-1, 1, nbpoints)
X, Y = np.meshgrid(c,c)
points = np.column_stack((X.ravel(), Y.ravel()))
for p in range(3):
    if (p==0): 
        titre =  "norme $\|.\|_1$"
        d = np.array([point for point in points if np.abs((np.linalg.norm(point, ord=1)-1))<0.01])
        x = d[:,0]
        y = d[:,1]
    elif (p==1): #norme 2
        titre = "norme $\|.\|_2$"
 #       d = np.array([point for point in points if np.linalg.norm(point, ord=2) == 1])

        theta = np.linspace(0, 2*np.pi, nbpoints)
        x = np.cos(theta)
        y = np.sin(theta)
    else:
        titre =  "norme $\|.\|_\infty$"
        d = np.array([point for point in points if np.linalg.norm(point, ord=np.inf) == 1])
        x = d[:,0]
        y = d[:,1]


    vec = np.array([x,y])
    A = np.array([[1,2],[0,2]])
    res = np.matmul(A,vec)
    if (p==0):
        vmax = [0,1]

    elif (p==1):
        l,v = LA.eig(np.matmul(A.transpose(),A))
        vmax = v[:,1]
        
    else:
        vmax = [1,1]
        
    Avmax = np.matmul(A,vmax)


    ax[p][0].spines['left'].set_position('center')
    ax[p][0].spines['bottom'].set_position('center')
    ax[p][0].spines['right'].set_color('none')
    ax[p][0].spines['top'].set_color('none')
    ax[p][0].xaxis.set_ticks(np.arange(-1, 2,1.0))
    ax[p][0].yaxis.set_ticks(np.arange(-1, 3,1.0))
    ax[p][0].plot([1],[0],marker='X',c='r',ms=10,label=())
    ax[p][0].plot([0],[1],marker='X',c='g',ms=10)
    ax[p][0].plot([vmax[0]],vmax[1],marker='*',c='b',ms=10)
    ax[p][0].plot([0,vmax[0]], [0,vmax[1]], linestyle='--', color='b')
    ax[p][0].scatter(x,y,s=1)
    ax[p][0].axis('equal')  

    ax[p][1].spines['left'].set_position('center')
    ax[p][1].spines['bottom'].set_position('center')
    ax[p][1].spines['right'].set_color('none')
    ax[p][1].spines['top'].set_color('none')
    ax[p][1].xaxis.set_ticks(np.arange(-2, 3,1.0))
    ax[p][1].yaxis.set_ticks(np.arange(-3, 4,1.0))
    ax[p][1].xaxis.set_ticks_position('bottom')
    ax[p][1].yaxis.set_ticks_position('left')
    ax[p][1].plot([1],[0],marker='X',c='r',ms=10)
    ax[p][1].plot([2],[2],marker='X',c='g',ms=10)
    ax[p][1].plot([Avmax[0]],Avmax[1],marker='*',c='b',ms=10)
    ax[p][1].plot([0,Avmax[0]], [0,Avmax[1]], linestyle='--', color='b')
    ax[p][1].axhline(0, color='black',linewidth=0.5)
    ax[p][1].axvline(0, color='black',linewidth=0.5)
    ax[p][1].scatter(res[0],res[1],s=1)
    ax[p][1].axis('equal')  
    ax[p][1].set_title('Effet de $A$ sur la boule unité, '+ titre)

plt.tight_layout()
plt.savefig('effetnorme.png',dpi=200)
plt.show()





