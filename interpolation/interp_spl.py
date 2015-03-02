import numpy as np

## Cette fonction permet de calculer les valeurs (ZO) aux points 
## d'interet par une interpolation par splines.

def interp_spl(X,Y,Zalti,Zinterp,p):

    tail=X.shape[0]                              ## nombre de lignes du vecteur X.
    X=X[:,np.newaxis]
    Y=Y[:,np.newaxis]
    Zalti=Zalti[:,np.newaxis]
    A=np.ones((tail,3))                          ## matrice formee des 3 vecteurs X, Y et V1.
    A[:,1]=X[:,0]
    A[:,2]=Y[:,0]
    B=np.zeros((tail,tail))                      ## initialisation d'une matrice de tail lignes et tail colonnes. 
    for i in range(0,tail):
        Bi = np.sqrt(np.power(X-X[i],2)+np.power(Y-Y[i],2)+np.power(Zalti-Zalti[i],2))
        Bi = np.squeeze(Bi,1)                     #removes single-dimensional entries form the shape
        for i in range(0,Bi.shape[0]):
            if Bi[i] == 0: Bi[i] = 1e-10
        B[:,i] = Bi * np.log(np.sqrt(Bi))        ## remplissage colonne par colonne de B.
        B[i,i] = p                               ## elements diagonaux de B sont nuls.

    C=np.zeros((3,3))                            ## initialisation de la matrice C de 3 lignes et 3 colonnes.
    D = A.transpose()
        
    M=np.bmat('A,B;C,D')                         ## creation de la matrice M qui est la concatenation des matrices A, B , C et A'.
    Zinterp=np.append(Zinterp,[0,0,0])
    Zinterp=Zinterp[np.newaxis,:]
    nul=np.zeros((1,Zinterp.shape[1]))
    Z1=np.bmat('Zinterp;nul;nul;nul').transpose()      ## creation du vecteur Z1 contenant les valeurs Z et dont les 3 dernieres lignes sont nulles.

    SORT=np.linalg.pinv(M.transpose().dot(M)).dot(M.transpose()).dot(Z1)                  ##obtention du vecteur des parametres (a0,a1,a2,b1,...,bn) par moindres carres.
    SORT=np.squeeze([x for x in SORT[:,0]])

    '''
    d = np.hypot(X-X0[i,k],Y-Y0[i,k])   ## calcul des distances entre les points de coordonnees (X,Y) et les points de coordonees (X0,Y0).
    K = np.squeeze(d*np.log(np.sqrt(d)))            ## calcul de la fonction K
    Z0[i,k] = SORT[0] + SORT[1]*X0[i,k] + SORT[2]*Y0[i,k] + sum(SORT[3:]*K)  ## calcul d'une ligne du vecteur Z0.
    '''
    
    return SORT
