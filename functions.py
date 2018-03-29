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
    "time measure function"
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

################
# Genomes size #
################

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
    """formatting numeric values
    input : string with numeric value and a comma as thousand separator
    output : int or float"""
    try:
        return int(a.replace(",",""))
    except:
        return float(a.replace(",",""))

def motif(bact="Pa14"): #position - direction - name - mismatch
    """Reading in an associated file patterns position 
    input : file name (in "DATA/YRS_%s.csv")
    output : list of tuples with:
        - end of patterns position
        - bool with patterns orientation
        - mismatches count"""
    doc=open("DATA/YRS_%s.csv"%bact,"rb")
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


def pourcentMotif(insert,bact="Pa14"):
    """rate of the inserts close to patterns (define with motif())
    input : the inserts and bacteria
    output : float"""
    YRS=motif(bact)
        
    z=0.
    for yrs in YRS:
        i,sens,mismatch=yrs
        if sens: 
            z+=np.sum(insert[1:,i+49:i+81])
        elif not sens:
            z+=np.sum(insert[1:,i-81:i-49])
    a=z/np.sum(insert[1:,:])
    #return round(a*100,3)
    return a

def NsiI():
    """"Reading in an associated file NsiI restriction sites positions (in DATA/NsiI.csv)
    output : list of tuples with:
        - beginning of NsiI position
        - end of NsiI position"""
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
    """from NsiI() made list of all NsiI restriction sites positions (and one on each side)
    output : list of int"""
    N=[]
    for i in NsiI():
        N+=range(i[0]-1,i[1]+2)
    return N

def genePAO1():
    """Reading in an associated file genes positions for PAO1 genome =
    output : list of tuples with:
        - beginning of gene
        - end of gene"""
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
    """from genePAO1() made list of all genes positions
    input : cutoff
    output : list of int"""
    N=[]
    for i in genePAO1():
        N+=range(i[0]+cut,i[1]+1-cut)
    return list(set(N))

def inter_list(cut=0):
    """made list of all intergenic positions
    input : cutoffs
    output : liste of int"""
    return list(set(range(1,Size["PAOI"]))-set(gene_list(cut=cut)))

def addition(liste):
    """sum of 2 or 3 inserts list
    input : lists
    output : summoned list"""
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
    """inserts clustering by a variable binding step
    input : inserts list and binding step
    output : bind inserts list"""
    g=np.zeros((3,int(ceil(float(len(insert[0]))/pas))),dtype=int)
    
    for i in range(len(g[0])):
        g[:,i]=[insert[0,pas*i],sum(insert[1,pas*i:pas*(i+1)]),sum(insert[2,pas*i:pas*(i+1)])]
    return g

def selectinrange(dico,r):
    """select inserts in a range
    input : inserts dict and positions to keep
    output : inserts dict"""
    out={}
    
    for k in dico.keys():
        if dico[k][0] in r:
            out[k]=dico[k]
    return out

def readsName(string):
    """cleanning reads name
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
    """Reading and formatting inserts in an attached file (product with Geneious)
    input : path to csv file 
    output : a tuple with 2 dictionaries (Beginning inserts and End inserts). Each dictionary keeps read name for key and tuple (position, orientation) for values"""
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
    """From inserts dict, count inserts present along the genome for each direction
    input : inserts dict and genome size
    output : tree lines numpy array:
        - position
        - forward inserts count
        - reverse inserts count"""
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
    """Analysed end of encapsidation
    input : 
        - dict of beginning insertions
        - dict of end insertions
        - genome size
        - bool for exclusion  of reads close to NsiI_list()
    output : tuple of 3 results
        - list of all encapsidated fragments size
        - end encapsidation inserts
        - dict of unused reads"""
    if Nsi:
        nsiI=NsiI_list()
    else:
        nsiI=[]
    
    
    Reads=0
    Taille=[]
    cut=0
    
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

###################################################
# Functions - graph() and other graphic functions #
###################################################

def graph(insert,name="graph",ticks=None,titre=None,NSI=None,save=False,size=None,Ymax=None):
    """ main graphical function
    input : 
        - vector of reads counts
        - graph name
        - ticks list
        - graph title
        - list of NsiI site to plot
        - bool for saving (or plotting)
        - graph size
        - Ymax
    output : graph"""
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
    """mirror graph with forward inserts (up) and reverse inserts (bottom)
    input : like graph()
    output : graph"""
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

def graph_multiline(insert,line=10,ticks=200,name="graph",titre=None,NSI=None,save=False,size=(16,2.8)):
    """display a multiline graph
    input : like graph()
        - line : number of lines to display
    output : graphs"""
    pas=ceil(float(len(insert[0])/line))
    Ymax=insert[1:,:].sum(axis=0).max()
    for i in range(line):        
        graph(insert[:,i*pas:(i+1)*pas],ticks=ticks,name="%s_%ion%i"%(name,i+1,line),titre=titre,NSI=NSI,save=save,size=size,Ymax=Ymax)

def graph_pics(insert,pic,delta=25,ticks=25,name="pic",titre=None,NSI=None,save=False,size=None):
    """display graphs each time sums of inserts exceeds a variable threshold
    input : like graph()
        - pic : threshold
        - delta : range around peaks
    output : graphs"""
    insertSum=insert[1:,:].sum(axis=0)
    big=insertSum>pic
    for i in range(len(big)):
        if big[i]:
            graph(insert[:,i-delta:i+delta],ticks=ticks,name="%s_%i"%(name,i),titre=titre,NSI=NSI,save=save,size=size)

def TablePics(insert,pic=1000,delta=25,bact="Pa14"):
    """Build a CSV table and a GFF file focus on peaks
    input : 
        - vector of reads counts
        - peaks threshold
        - base around peaks
        - bacteria utilized
    output : CSV table GFF file"""
    insertSum=insert[1:,:].sum(axis=0)
    big=insertSum>pic
    
    YRS=np.array(motif(bact),dtype=int)
    NSI=np.array(NsiI(),dtype=int)
    
    fichier=open("fig/tableHotspot.csv", "wb")
    writer=csv.writer(fichier)
    writer.writerow(("Position","Nb Inserts","Distance Motif","Position Motif","Distance NsiI","Ratio forward"))
    
    GFF=open("fig/hotspot.gff", "wb")
    
    for i in range(len(big)):
        if big[i]:
            line=[]
            line.append(i)
            line.append(insertSum[i])
            
            closeMotif=YRS[abs(i-YRS[:,0])==min(abs(i-YRS[:,0])),0]
            line.append(abs(i-closeMotif))
            line.append(closeMotif)
            
            nsi1=abs(i-max(filter(lambda x:i-x>0,NSI[:,1])))
            nsi2=abs(i-min(filter(lambda x:i-x<0,NSI[:,0])))
            line.append((nsi1,nsi2))
            
            line.append(float(np.sum(insert[1,i-delta:i+delta]))/sum(insertSum[i-delta:i+delta]))
            writer.writerow(line)
            
            GFF.write("%s\tPython\thotspot\t%i\t%i\t.\t+\t.\t\n"%(bact,i,i))
            
                        
    fichier.close()
    GFF.close()

#################
# others graphs #
#################

def graph_separate(insert,save=False):
    """display two separeted graphs : one for forward inserts and one for reverse
    input : vector of reads counts
    output : graphs"""
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
    """display rapidly a global graph by adding forward and reverse inserts
    input : vector of reads counts, graph name, graph title, bool for saving
    output : graph"""
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

############################
# Histogram of insert size #
############################

def histo(taille,bins=None,name="histo",axes=None,titre=None,save=False,size=(16,8),Xmax=None):
    """histogram of insert size and binning if desired
    input : 
        - list of insert size (product with dataAnalysePaired())
        - list of binning positions
        - axes = [Xmin, Xmax, Ymin, Ymax]
        - graph title
        - bool for saving (or plotting)
        - graph size
        - ordinate Xmax
    output : histogram"""
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
        
    if axes:
        plt.axis(axes)
    
    #Save ?
    if save:
        plt.savefig("fig/%s.jpeg"%name,dpi=500,transparant=True)
    else:
        plt.show()
        
    plt.close()
    return

def histoDouble(taille1,taille2,bins1=None,bins2=None,name="histo",axes=None,titre=None,save=False,size=(16,8),Xmax=None):
    """double histogram
    input : like histo () with taille1 et taille2 for each inserts size and bin1 et bin2 associated
    output : histogram"""
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
        
    if axes:
        plt.axis(axes)

    
    #Save ?
    if save:
        plt.savefig("fig/%s.jpeg"%name,dpi=500,transparant=True)
    else:
        plt.show()
        
    plt.close()
    return
