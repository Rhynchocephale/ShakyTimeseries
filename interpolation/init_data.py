import numpy as np

## Chargement des donnees d'observation a interpoler, du trait de cote et creation de la grille vide
##---ENTREES---
##pasGrille : pas de la grille a creer
##dim: 0 for corrections in e, 1 for n, 2 for u
##---SORTIES---
##X : vecteur contenant les abscisses des observations
##Y : vecteur contenant les ordonnees des observations
##Z : vecteur contenant les observations
##trait_cote : matrice contenant les points du trait de cote
##X0 : matrice contenant les abscisses de la grille (xmin et xmax determines a partir du lot de donnees, pas passe en parametre)
##Y0 : matrice contenant les ordonnees de la grille (ymin et ymax determines a partir du lot de donnees, pas passe en parametre)
def init_data(enu,pos,step,dim,day):
    
    X=pos[:,0]
    Y=pos[:,1]
    Z=enu[:,day,dim]

    # Creation de la grille vide
    [X0,Y0] = np.meshgrid(range(int(np.floor(min(X))),int(np.floor(max(X))+1),step),range(int(np.floor(min(Y))),int(np.floor(max(Y))+1),step))
    
    np.savetxt('X',X)
    np.savetxt('Y',Y)
    np.savetxt('Z',Z)
    np.savetxt('X0',X0)
    np.savetxt('Y0',Y0)
    
    return [X,Y,Z,X0,Y0]
    
def init_data_nogrid(enu,pos,dim,day):

    is_bad=enu[:,day,3]
    zero_indices = np.nonzero(is_bad == 0)[0]
    nonzero_indices = np.nonzero(is_bad)[0]
    
    Zinterp=enu[:,day,dim]    
    X=np.array([x[0] for x in pos.values()])
    Y=np.array([x[1] for x in pos.values()])
    Zalti=np.array([x[2] for x in pos.values()])

    
    X=X[zero_indices]
    Y=Y[zero_indices]
    Zalti=Zalti[zero_indices]
    Zinterp=Zinterp[zero_indices]
        
    return [X,Y,Zalti,Zinterp,nonzero_indices]

