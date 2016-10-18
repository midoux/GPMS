#! /usr/bin/env python
# -*- coding: utf-8 -*-

##########
# Import #
##########

import pdb

import os
import numpy as np
import csv
import matplotlib.pyplot as plt
import time
from math import ceil

csv.field_size_limit(10000000)

def t():
    "fonction pour mesure du temps de calcul"
    try:
        t1=time.time()-t0
        if t1<60:
            print "Time : %0.2f s."%(t1)
        else:
            print "Time : %i:%2.0f"%(int(t1)/60, t1%60)
        return
    except:
        return
    
t0=time.time()

#####################################
#Definition des tailles des génomes #
#####################################

Size={}
Size["PAOI"]=6264404 #PAOIw
Size["Pa14"]=6541482 #PA14or
Size["II10"]=6288645
Size["Ab30"]=37205
Size["2P1"]=37087
Size["PM105"]=39593

#############
# Functions #
#############

def numeric(a):
    """mise en forme des valeurs numériques
    input : string comprenant une valeur numérique et des virgules comme séparateur de millier
    output : int"""
    try:
        return int(a.replace(",",""))
    except:
        return float(a.replace(",",""))

def motif(bacterie="Pa14"): #position - direction - name - mismatch
    """Lecture dans le fichier associé des positions des motifs
    input : nom du fichier (positionné dans "DATA/YRS_%s.csv")
    output : liste avec des tuples comprenant :
        - la position de fin du motif
        - un bool de la direction du motif
        - le nombre de mismatch"""
    doc=open("DATA/YRS_%s.csv"%bacterie,"rb")
    table=csv.reader(doc)
    
    intitule=table.next()
    Nam=intitule.index("Name")
    Min=intitule.index("Minimum")
    Max=intitule.index("Maximum")
    Dir=intitule.index("Direction")
    Mis=intitule.index("Mismatches")
    
    X=[]
    
    try:
        while 1:
            data=table.next()
            x=[]
            
            if data[Dir]=='forward':
                x.append(numeric(data[Max]))
                x.append(True)
            elif data[Dir]=='reverse':
                x.append(numeric(data[Min]))
                x.append(False)
            else:
                print data
                
            x.append(data[Mis])
                
            X.append(x)
            
    except StopIteration:
        pass
    
    doc.close()
    
    return X


def pourcentMotif(insert,bacterie="Pa14"):
    """ratio d'insert présents à proximité des motifs définies avec la fonction motif()
    input : les inserts et le noms de la bactérie associé
    output : float du ratio"""
    YRS=motif(bacterie)
        
    z=0.
    for yrs in YRS:
        i,sens,mismatch=yrs
        if sens: 
            #z+=np.sum(insert[1:,i+52:i+78])
            z+=np.sum(insert[1:,i+49:i+81])
        elif not sens:
            #z+=np.sum(insert[1:,i-78:i-52])
            z+=np.sum(insert[1:,i-81:i-49])
    a=z/np.sum(insert[1:,:])
    #return round(a*100,3)
    return a

def NsiI():
    """Lecture dans le fichier associé des positions des sites de restriction NsiI (positionné dans DATA/NsiI.csv)
    output : Liste de tuples comprenant le début et la fin de chaque site NsiI"""
    doc=open("DATA/NsiI.csv","rb")
    table=csv.reader(doc)
    
    intitule=table.next()
    Nam=intitule.index("Name")
    Min=intitule.index("Minimum")
    Max=intitule.index("Maximum")
    
    X=[]
    
    try:
        while 1:
            data=table.next()
            X.append((numeric(data[Min]),numeric(data[Max]))) 
    
    except StopIteration:
        pass
    
    doc.close()
    
    return X

def NsiI_list():
    """ a partir de la liste de la fonction NsiI() produit un liste comprenant l'ensemble des positions des sites NsiI (et une base de part et d'autres)
    output : liste de int des positions"""
    N=[]
    for i in NsiI():
        N+=range(i[0]-1,i[1]+2)
    return N

def genePAO1():
    """Lecture des positions des genes a partir du fichier "DATA/gene_PAOI.csv "
    output : liste des tuples des int des positions des gènes"""
    doc=open("DATA/gene_PAOI.csv","rb")
    table=csv.reader(doc)
    
    intitule=table.next()
    Min=intitule.index("Minimum")
    Max=intitule.index("Maximum")
    
    X=[]
    
    try:
        while 1:
            data=table.next()
            X.append((numeric(data[Min]),numeric(data[Max]))) 
    
    except StopIteration:
        pass
    
    doc.close()
    
    return X

def gene_list(cut=0):
    """Génère une liste de toute les positions intragéniques a partir de la fonction précédente
    input : cutoff des positions
    output : liste de int des positions"""
    N=[]
    for i in genePAO1():
        N+=range(i[0]+cut,i[1]+1-cut)
    return list(set(N))

def inter_list(cut=0):
    """Génère une liste de toute les positions intergénique
    input : cutoff des positions
    output : liste de int des positions"""
    return list(set(range(1,Size["PAOI"]))-set(gene_list(cut=cut)))

def addition(liste):
    """somme une liste (2 ou 3) d'inserts
    input : liste
    output : inserts sommés"""
    if len(liste)==2:
        if list(liste[0][0])==list(liste[1][0]):
            return np.array([liste[0][0],liste[0][1]+liste[1][1],liste[0][2]+liste[1][2]])
    elif len(liste)==3:
        if list(liste[0][0])==list(liste[1][0]):
            if list(liste[1][0])==list(liste[2][0]):
                return np.array([liste[0][0],liste[0][1]+liste[1][1]+liste[2][1],liste[0][2]+liste[1][2]+liste[2][2]])
    else:
        print "ERR0R" 

def group(insert,pas=1):
    """regroupement des inserts suivant un pas variable
    input : inserts et pas de binding
    output : inserts regroupés"""
    g=np.zeros((3,ceil(float(len(insert[0]))/pas)),dtype=int)
    
    for i in range(len(g[0])):
        g[:,i]=[insert[0,pas*i],sum(insert[1,pas*i:pas*(i+1)]),sum(insert[2,pas*i:pas*(i+1)])]
    return g

def selectinrange(dico,r):
    """sélection uniquement des inserts dans une liste de positions données
    input : dictionnaire des inserts et liste des positions à conserver
    output : dictionnaire des inserts conserver"""
    out={}
    
    for k in dico.keys():
        if dico[k][0] in r:
            out[k]=dico[k]
    return out

def readsName(string):
    """parseur des noms de reads
    input : string
    output : string"""
    if "_" in string:
        return string.split("_")[0]
    elif "/" in string:
        return string.split("/")[0]
    elif " " in string:
        return string.split(" ")[0]
    else:
        return string

def csvData (document):
    """Lecture et mise en forme des inserts présent dans un fichier joint (produit avec Geneious)
    input : nom du fichier csv a lire 
    output : un tuple contenant 2 dictionnaires (Debut et Fin d'insertion). Chaque dictionnaire prend comme clé le nom du read parsé et comme valeur un tuple avec un int pour la position et un bool pour l'orientation"""
    #Dico avec nom du read en clé et tuples (int_position , bool_direction)
    Start={}
    Stop={}
    
    D=open(document,"rb")
    table=csv.reader(D)
    
    intitule=table.next()
    Typ=intitule.index("Type")
    Min=intitule.index("Minimum")
    Dir=intitule.index("Direction")
    From=intitule.index("Transferred From")
    
    
    try:
        while 1:
            data=table.next()
            
            if data[Typ]=="start":
                for read in data[From].split(", "):
                    if Start.has_key(readsName(read)):
                        if Start[readsName(read)]!=(numeric(data[Min]), data[Dir]=='forward'):
                            if data[Dir]==Start[readsName(read)][1]:
                                if data[Dir]:
                                    Start[readsName(read)]=(min(Start[readsName(read)][0],numeric(data[Min])), data[Dir]=='forward')
                                else :
                                    Start[readsName(read)]=(max(Start[readsName(read)][0],numeric(data[Min])), data[Dir]=='forward')
                            else:
                                del Start[readsName(read)]
                    else:
                        Start[readsName(read)]=(numeric(data[Min]), data[Dir]=='forward')
            
            elif data[Typ]=="stop":
                for read in data[From].split(", "):
                    if Stop.has_key(readsName(read)):
                        Stop[readsName(read)].append(((numeric(data[Min]), data[Dir]=='forward')))
                    else:
                        Stop[readsName(read)]=[(numeric(data[Min]), data[Dir]=='forward')]
                        
            else:
                print data
    
    except StopIteration:
        pass
    
    D.close()
    
    return (Start, Stop)

def dataAnalyseStart (Start,genomeSize):
    """Décompte des insertions présentes dans chaque sens le long du génome en fonction du dictionnaire des inserts
    input : dictionnaire des insert et taille du génome
    output : array de int de trois lignes :
        - indice de la position
        - nombre d'insertions forward
        - nombre d'insertions reverse"""
    insert=np.array([np.arange(1,genomeSize+1,dtype=int),
           np.zeros(genomeSize,dtype=int),
           np.zeros(genomeSize,dtype=int)])
    
    for i in Start:
        if Start[i][1]:
            insert[1,Start[i][0]-1]+=1
        else:
            insert[2,Start[i][0]-1]+=1
    
    return insert

def dataAnalysePaired (Start,Stop,genomeSize,Nsi=False):
    """Analyse des jeux d'insertion avec sites de fin d'encapsidation
    input : 
        - dictionnaire des débuts d'insertions
        - dictionnaire des fins d'insertions
        - taille du genome
        - bool pour l'exclusion des fins d'insertions a proximité d'un site de restriction définie dans NsiI_list()
    output : tuple de 3 résultats
        - liste d'int de toutes les taille des fragments encapsidés
        - insert de fin d'encapsidation
        - dictionnaire des reads non utilisés"""
    if Nsi:
        nsiI=NsiI_list()
    else:
        nsiI=[]
    
    
    Reads=0
    Taille=[]
    cut=0 #Nb de reads coupés par NsiI
    
    unused={}
    
    insert=np.array([np.arange(1,genomeSize+1,dtype=int),
           np.zeros(genomeSize,dtype=int),
           np.zeros(genomeSize,dtype=int)])
    
    for i in Start.keys():
        if Stop.has_key(i):
            Reads+=1
            (sPos,sDir)=Start[i]
            
            end=Stop[i]
            
            if type(end)==list:
                if len(end)==1:
                    if end[0][1]!=sDir:
                        ePos=end[0][0]
                    elif end[0][1]==sDir:
                        ePos=None
                    
                elif len(end)>1:
                    for e in range(len(end)):
                        if end[e][1]!=sDir:
                            end[e]=end[e][0]
                        elif end[e][1]==sDir:
                            end[e]=None
                
                    if sDir:
                        ePos=max(end)
                    else:
                        ePos=min(x for x in end if x is not None)
                        
                    Stop[i]=(ePos,-sDir)
            elif type(end)==tuple:
                (ePos,eDir)=Stop[i]
                        
            if ePos in nsiI:
                cut+=1
                
            else:
                if ePos:
                    taille=min(abs(sPos-ePos)+1,genomeSize-abs(sPos-ePos)+1)
                    if taille<3000:
                        Taille.append(taille)
                        if sDir:
                            insert[1,ePos-1]+=1
                        elif not sDir:
                            insert[2,ePos-1]+=1
                            
                    else:
                        unused[i]=Stop[i]
        
    print "Reads utilisés = %i"%Reads
    if Nsi:
        print "Reads lysés par NsiI = %0.2f %%"%(100*(float(cut)/Reads))
    print "Reads Fail = %0.2f %%"%(float(100*(len(unused))/Reads))
    
    return (Taille, insert, unused)

#############################################
# Functions - graph() et graph en découlant #
#############################################

def graph(insert,name="graph",ticks=None,titre=None,NSI=None,save=False,size=None,Ymax=None):
    """ fonction principale pour les graph d'insertion
    input : 
        - vecteur du decompte des reads
        - nom du graph (pour sauvegarde)
        - liste des ticks a affiher
        - titre du graph
        - liste des site NsiI a afficher
        - bool de sauvegarde (sinon affichage)
        - taille du graph
        - Ymax
    output : sauvegarde ou affichage du graph"""
    if not titre:
        titre=name
        
    #plt.clf()
    fig = plt.figure(figsize=size)
    
    p1 = plt.bar(insert[0,:],
                 insert[1,:],
                 width=np.min(insert[0,1:]-insert[0,:-1]),
                 color='r',
                 edgecolor=(0, 0, 0, 0.0))
    
    p2 = plt.bar(insert[0,:],
                 insert[2,:],
                 width=np.min(insert[0,1:]-insert[0,:-1]),
                 color='b',
                 edgecolor=(0, 0, 0, 0.0),
                 bottom=insert[1,:])
    
    
    plt.title(titre)
    plt.legend( (p1[0], p2[0]), ('forward', 'reverse') )
    plt.xlabel("Bacterial position")
    plt.ylabel("Inserts count")
    
    ind=insert[0]
    if ticks:
        plt.xticks(ind[ind%ticks==0],ind[ind%ticks==0])
    plt.xlim(insert[0,0],insert[0,-1])
    if Ymax:
        plt.ylim(0,Ymax)
    else:
        Ymax=insert[1:,:].sum(axis=0).max()
    
    #Site NsiI
    if NSI:
        for nsi in NSI:
            markerline, stemlines, baseline = plt.stem(np.array([nsi]), np.array([Ymax]),':')
    #Save ?
    if save:
        plt.savefig("fig/%s.jpeg"%name,dpi=600,transparant=True)
    else:
        plt.show()
        
    plt.close()
    return

def graphMirror(insert,name="graph",ticks=None,titre=None,NSI=None,save=False,size=None,Ymax=None):
    """ affichage du graph en mirroir avec les inserts forwards vers le haut et les reverses vers le bas
    input : identique à graph()
    output : sauvegarde ou affichage du graph"""
    if not titre:
        titre=name
        
    #plt.clf()
    fig = plt.figure(figsize=size)
    
    p1 = plt.bar(insert[0,:],
                 insert[1,:],
                 width=np.min(insert[0,1:]-insert[0,:-1]),
                 color='r',
                 edgecolor=(0, 0, 0, 0.0))
    
    p2 = plt.bar(insert[0,:]-4,
                 -insert[2,:],
                 width=np.min(insert[0,1:]-insert[0,:-1]),
                 color='b',
                 edgecolor=(0, 0, 0, 0.0))
    
    
    plt.title(titre)
    plt.legend( (p1[0], p2[0]), ('forward', 'reverse') )
    plt.xlabel("Bacterial position")
    plt.ylabel("Inserts count")
    
    ind=insert[0]
    if ticks:
        plt.xticks(ind[ind%ticks==0],ind[ind%ticks==0])
    plt.xlim(insert[0,0],insert[0,-1])
    if Ymax:
        plt.ylim(-Ymax,Ymax)
        YmaxUp=Ymax
        YmaxDown=-Ymax
    else:
        YmaxUp=insert[1,:].max()
        YmaxDown=-insert[2,:].max()
    
    #Site NsiI
    if NSI:
        for nsi in NSI:
            markerline, stemlines, baseline = plt.stem(np.array([nsi]), np.array([YmaxUp]),':')
            markerline, stemlines, baseline = plt.stem(np.array([nsi]), np.array([YmaxDown]),':')
    #Save ?
    if save:
        plt.savefig("fig/%s.jpeg"%name,dpi=600,transparant=True)
    else:
        plt.show()
        
    plt.close()
    return

def graph_multiline(insert,ligne=10,ticks=200,name="graph",titre=None,NSI=None,save=False,size=(16,2.8)):
    """ affichage du graph découpé en plusieurs ligne
    input : identique à graph()
        - ligne : nombre de lignes à afficher
    output : sauvegarde ou affichage des graphs"""
    pas=ceil(float(len(insert[0])/ligne))
    Ymax=insert[1:,:].sum(axis=0).max()
    for i in range(ligne):        
        graph(insert[:,i*pas:(i+1)*pas],ticks=ticks,name="%s_%ion%i"%(name,i+1,ligne),titre=titre,NSI=NSI,save=save,size=size,Ymax=Ymax)

def graph_pics(insert,pic,delta=25,ticks=25,name="pic",titre=None,NSI=None,save=False,size=None):
    """ affichage des graphs au niveau des positions ou le nombre d'inserts dépasse un certain seuil
    input : identique à graph()
        - pic : seuil a partir duquel les graphs sont affichés
        - delta : nombre de base de part et d'autres des pics
    output : sauvegarde ou affichage du graph"""
    insertSum=insert[1:,:].sum(axis=0)
    big=insertSum>pic
    for i in range(len(big)):
        if big[i]:
            graph(insert[:,i-delta:i+delta],ticks=ticks,name="%s_%i"%(name,i),titre=titre,NSI=NSI,save=save,size=size)

def TablePics(insert,pic=1000,delta=25,bacterie="Pa14"):
    """ Generation d'une table CSV et d'une annotation GFF des pics d'insert
    input : 
        - vecteur du decompte des reads
        - seuil a partir duquel les graphs sont affichés
        - nombre de base de part et d'autres des pics
        - bacterie utilisée
    output : enregistrement de la table CSV et du fichier GFF"""
    insertSum=insert[1:,:].sum(axis=0)
    big=insertSum>pic
    
    YRS=np.array(motif(bacterie),dtype=int)
    NSI=np.array(NsiI(),dtype=int)
    
    fichier=open("fig/tableHotspot.csv", "wb")
    writer=csv.writer(fichier)
    writer.writerow(("Position","Nb Inserts","Distance Motif","Position Motif","Distance NsiI","Ratio forward"))
    
    GFF=open("fig/hotspot.gff", "wb")
    
    for i in range(len(big)):
        if big[i]:
            ligne=[]
            ligne.append(i) #Position
            ligne.append(insertSum[i]) #Nb Inserts
            
            closeMotif=YRS[abs(i-YRS[:,0])==min(abs(i-YRS[:,0])),0]
            ligne.append(abs(i-closeMotif)) #Distance Motif
            ligne.append(closeMotif)#Position Motif
            
            nsi1=abs(i-max(filter(lambda x:i-x>0,NSI[:,1])))
            nsi2=abs(i-min(filter(lambda x:i-x<0,NSI[:,0])))
            ligne.append((nsi1,nsi2)) #Distance NsiI
            
            ligne.append(float(np.sum(insert[1,i-delta:i+delta]))/sum(insertSum[i-delta:i+delta])) #Ratio forward
            
            writer.writerow(ligne)
            
            GFF.write("%s\tPython\thotspot\t%i\t%i\t.\t+\t.\t\n"%(bacterie,i,i))
            
                        
    fichier.close()
    GFF.close()

##############################
# Functions - autres graphes #
##############################

def graph_separate(insert,save=False):
    """ affichage de 2 graphs, un pour les inserts forwards et l'autre pour les inserts reverses
    input : vecteur du decompte des reads et bool de sauvegarde
    output : sauvegarde ou affichage des graphs"""
    forward=insert[1,:]
    reverse=insert[2,:]
    
    plt.clf()
    m=max([max(reverse),max(forward)])
    M=0
    while M<m:
        M+=100

    plt.figure(1)
    plt.xlim(0,len(forward))
    plt.ylim(0,M)
    plt.title('Inser Forward')
    plt.plot(forward,color='r')
    plt.xlabel("Bacterial position")
    plt.ylabel("Inserts count")
    if save:
        plt.savefig("fig/graph_Forward.jpeg",dpi=500,transparant=True)
    else:
        plt.show()
    
    plt.figure(2)
    plt.xlim(0,len(reverse))
    plt.ylim(0,M)
    plt.title('Inser Reverse')
    plt.plot(reverse,color='b')
    plt.xlabel("Bacterial position")
    plt.ylabel("Inserts count")
    if save:
        plt.savefig("fig/graph_Reverse.jpeg",dpi=500,transparant=True)
    else:
        plt.show()
        
    plt.close()
    return

def graph_global(insert,name="graph",titre=None,save=False):
    """affichage d'un graph en sommant les inserts forwards et reverse
    input : vecteur du decompte des reads, nom du graph, titre du graph et bool de sauvegarde
    output : sauvegarde ou affichage du graph"""
    if not titre:
        titre=name
        
    plt.clf()
    plt.title(titre)
    plt.plot(insert[1:,:].sum(axis=0))
    plt.xlim(0,len(insert[0]))
    plt.xticks([0,1000000,2000000,3000000,4000000,5000000,6000000],["0M","1M","2M","3M","4M","5M","6M"])
    plt.xlabel("Bacterial position")
    plt.ylabel("Inserts count")
    if save:
        plt.savefig("fig/%s.jpeg"%name,dpi=500,transparant=True)
    else:
        plt.show()
        
    plt.close()
    return

#######################
# Histo taille insert #
#######################

def histo(taille,bins=None,name="histo",axies=None,titre=None,save=False,size=(16,8),Xmax=None):
    """histogramme de la taille des inserts et binning si souhaité
    input : 
        - liste de la taille des inserts produit avec la fonction dataAnalysePaired()
        - liste des position de binning
        - liste des axes [Xmin, Xmax, Ymin, Ymax]
        - titre du graph
        - bool de sauvegarde
        - taille du graph
        - ordonnée Xmax
    output : sauvegarde ou affichage de l'histogramme"""
    if not titre:
        titre=name
        
    plt.clf()
    
    fig, ax1 = plt.subplots(figsize=size)
    ax1.hist(taille,bins=max(taille)-min(taille),color="blue", histtype="step")
    
    ax1.set_xlabel('Inserts size')
    ax1.set_ylabel('Global distribution', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    
    if bins:
        ax2 = ax1.twinx()
        ax2.hist(taille,bins=bins,color="red", histtype="step")
        ax2.set_ylabel('Clustered distribution', color='r')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
    
    plt.title(titre)
    if Xmax:
        plt.xlim((0,Xmax))
        
    if axies:
        plt.axis(axies)

    
    #Save ?
    if save:
        plt.savefig("fig/%s.jpeg"%name,dpi=500,transparant=True)
    else:
        plt.show()
        
    plt.close()
    return

def histoDouble(taille1,taille2,bins1=None,bins2=None,name="histo",axies=None,titre=None,save=False,size=(16,8),Xmax=None):
    """histogramme de la taille de deux jeux d'inserts et binning si souhaité
    input : identique a histo () avec taille1 et taille 2 pour chaque jeux de taille d'inserts et bin1 et bin2 associé
    output : sauvegarde ou affichage de l'histogramme"""
    if not titre:
        titre=name
        
    plt.clf()
    
    fig, ax1 = plt.subplots(figsize=size)
    if bins1:
        ax1.hist(taille1,bins=bins1,color="blue", histtype="step")
    else:
        ax1.hist(taille1,bins=max(taille1)-min(taille1),color="blue", histtype="step")
    
    ax1.set_xlabel('Inserts size')
    ax1.set_ylabel('Global distribution', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    
    
    ax2 = ax1.twinx()
    if bins2:
        ax2.hist(taille2,bins=bins2,color="red", histtype="step")
    else:
        ax2.hist(taille2,bins=max(taille2)-min(taille2),color="red", histtype="step")
    ax2.set_ylabel('Focus distribution', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    
    plt.title(titre)
    if Xmax:
        plt.xlim((0,Xmax))
        
    if axies:
        plt.axis(axies)

    
    #Save ?
    if save:
        plt.savefig("fig/%s.jpeg"%name,dpi=500,transparant=True)
    else:
        plt.show()
        
    plt.close()
    return


