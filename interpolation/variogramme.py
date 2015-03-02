import numpy as np

def variogramme(X,Y,Zalti,Zinterp):
    
    lenx = len(X)
    p=(lenx-2)*(lenx-1)/2
    Tab=np.zeros(p,2)                  # matrice de p lignes et 2 colonnes
    k=0
    for i in range(0,lenx):
        for j in range(i+1,lenx):
            Tab[k,0]=np.sqrt(np.pow(X[i]-X[j],2)+np.pow(Y[i]-Y[j],2)+np.pow(Zalti[i]-Zalti[j],2))  # distance entre les points
            Tab[k,1]=(Zinterp[i]-Zinterp[j])^2                         # ecarts aux carres des variables regionalisees
            k+=1

    
    #figure(1)
    #plot (Tab(:,1),Tab(:,2),'+')
    '''IL SE PASSE DES TRUCS, ICI...'''
    maxi=max(Zinterp)     
    pas=max(Zinterp)-min(Zinterp)
    vario=zeros(maxi/pas,2)                                 # initialisation d'une matrice de max lignes et 2 colonnes.
    for i in range(0,maxi/pas):
    
        id=where(Tab[:,0]>=i*pas and Tab[:,1]<(i+1)*pas)     # cherche les indices des lignes de la premiere colonne de Tab ou les distances sont comprises entre (i-1)*pas et i*pas (par palier).
        
        vario[i,0]=mean(Tab[id,0])                         # moyenne des distances sur la premiere colonne de vario
        vario[i,1]=mean(Tab[id,1])/2                       # demi-moyenne des ecarts aux carres sur la deuxieme colonne de vario (estimateur de l'esperance)

    #gam0 calcule avec a0 et C0
    #gam=vario(:,2)
    #h=vario(:,1)
    #gam(h)=gam0(h) +(1-exp(-h^2/a0^2))*dC+da*C0*2*r^2*exp(-h^2/a0^2))/a0^3
    #O=A*P <=> Y=AX pour les moindres carres
    #P=[dC;da];
    #------------------Utilisation du modele exponentiel----------------------------
    a0=600000;
    C0=0.0035;                     # initialisation des parametres
    dC=0.0001; 
    da=40000;
    
    while (abs(dC/C0)>0.001 or abs(da/a0)>0.001):
        gam0=C0*(1-np.exp(-(vario[:,1]/a0)))                                           # on injecte les parametres initiaux dans gam0
        O=vario[:,1]-gam0
        vario0 = vario[:,0]
        A=[1-np.exp(-1*vario0/a0) -C0*vario0*np.exp(-pow(vario0/a0,2)/a0^2)/a0^2]     # remplissage de la matrice A=[dgam/dC,dgam/da] ligne par ligne des derivees partielles.
        P=np.linalg.pinv(A.transpose().dot(A)).dot(A.transpose()).dot(O)          # obtention du vecteur des parametres P=[dC;da]
        dC=P[0]
        da=P[1]
        a0=a0+da                                                                 # on incremente les parametres initiaux
        C0=C0+dC
    
    return [C0,a0,vario,gam0]
