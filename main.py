#! /usr/bin/env python
# -*- coding: utf-8 -*-


from functions import *

print "# Insertion de Ab30 dans un autre phage"

A_Ab30_DATA=csvData("DATA/A_Ab30_droite2.csv")
A_Ab30_insert=dataAnalyseStart(A_Ab30_DATA[0],Size["Ab30"])
A_Ab30_insert50=group(A_Ab30_insert,pas=50)

B_Ab30_DATA=csvData("DATA/B_Ab30_droite2.csv")
B_Ab30_insert=dataAnalyseStart(B_Ab30_DATA[0],Size["Ab30"])
B_Ab30_insert50=group(B_Ab30_insert,pas=50)


graph(A_Ab30_insert50,titre=u"A inserts distribution on Ab30 (group by 50bp)",name="f3_A",save=True)
graph(B_Ab30_insert50,titre=u"B inserts distribution on Ab30 (group by 50bp)",name="f3_B",save=True)

AB_Ab30_insert=addition([A_Ab30_insert,B_Ab30_insert])
graph(AB_Ab30_insert,titre=u"A+B inserts distribution on Ab30", name="f3_AB",save=True)

AB_Ab30_insert50=group(AB_Ab30_insert,pas=50)
graph(AB_Ab30_insert50,titre=u"A+B inserts distribution on Ab30 (group by 50bp)",name="f4",Ymax=400,save=True)

graph(AB_Ab30_insert[:,18499:19000],titre=u"A+B inserts distribution on Ab30 (focus on the mid)",name="f4_focus1",save=True)
graph(AB_Ab30_insert[:,9999:15000],titre=u"A+B inserts distribution on Ab30 (focus on 10k - 15k)",name="f4_focus2",save=True)
graph(AB_Ab30_insert[:,14999:20000],titre=u"A+B inserts distribution on Ab30 (focus on 15k - 20k)",name="f4_focus3",save=True)
graph(AB_Ab30_insert[:,17999:19000],titre=u"A+B inserts distribution on Ab30 (focus on 18k - 19k)",name="f4_focus4",save=True)

graph(AB_Ab30_insert[:,5399:7000],titre=u"A+B inserts distribution on Ab30 (focus on NsiI site)", name="f4_NsiI",NSI=range(6165,6170),save=True)

for banque in [A_Ab30_insert, B_Ab30_insert]:
    print "****"  
    print "inserts count =                    %i"%banque[1:,:].sum()
    print "ratio forward =                    %f"%(float(banque[1,:].sum())/banque[1:,:].sum())
    print "mean insertion =                   %f"%np.mean(banque[1:,:].sum(axis=0))
    print "variance =                         %f"%np.var(banque[1:,:].sum(axis=0))
print "****"

print "# Insertion of Ab30 in Pa14"

A_Pa14_DATA=csvData("DATA/A_Pa14_droite2.csv")
A_Pa14_insert=dataAnalyseStart(A_Pa14_DATA[0],Size["Pa14"])

graph_global(A_Pa14_insert,titre=u"A inserts distribution on PA14", name="f5_A",save=True)

B_Pa14_DATA=csvData("DATA/B_Pa14_droite2.csv")
B_Pa14_insert=dataAnalyseStart(B_Pa14_DATA[0],Size["Pa14"])

graph_global(B_Pa14_insert,titre=u"B inserts distribution on PA14", name="f5_B",save=True)

for banque in [A_Pa14_insert, B_Pa14_insert]:
    print "****"  
    print "inserts count =                    %i"%banque[1:,:].sum()
    print "ratio forward =                    %f"%(float(banque[1,:].sum())/banque[1:,:].sum())
    print "mean insertion =                   %f"%np.mean(banque[1:,:].sum(axis=0))
    print "variance =                         %f"%np.var(banque[1:,:].sum(axis=0))
    print "concentration closes to patterns = %f"%pourcentMotif(banque, "Pa14")
print "****"

AB_Pa14_insert=addition([A_Pa14_insert,B_Pa14_insert])

for banque in [AB_Pa14_insert]:
    print "****"  
    print "inserts count =                    %i"%banque[1:,:].sum()
    print "ratio forward =                    %f"%(float(banque[1,:].sum())/banque[1:,:].sum())
    print "mean insertion =                   %f"%np.mean(banque[1:,:].sum(axis=0))
    print "variance =                         %f"%np.var(banque[1:,:].sum(axis=0))
    print "concentration closes to patterns = %f"%pourcentMotif(banque, "Pa14")
print "****"

graph_global(AB_Pa14_insert,titre=u"A+B inserts distribution on Pa14", name="f6",save=True)

graph(AB_Pa14_insert[:,5217624:5217700],titre=u"Inserts distribution at a hotspot position",ticks=25,name="f7",NSI=NsiI_list(),save=True) 

graph_pics(AB_Pa14_insert,pic=2500,ticks=25,name="hotspot/",titre=None,NSI=NsiI_list(),save=True,size=None)

graph(AB_Pa14_insert[:,5729499:5736000],titre=u"Inserts distribution at two differents hotspot position",ticks=2000,name="f11a",NSI=NsiI_list(),save=True)
graph(AB_Pa14_insert[:,5729499:5736000],titre=u"Inserts distribution at two differents hotspot position",ticks=2000,name="f11b",save=True)

test=group(AB_Pa14_insert[:,5729499:5736000],pas=10)
graph(test,titre=u"Inserts distribution at two differents hotspot position (by 10bp)",ticks=2000,name="f11c",NSI=NsiI_list(),save=True)

print "# Excluding Hotspots"

AB_Pa14_insertBeg50=group(AB_Pa14_insert[:,0:40000],pas=50)
graph(AB_Pa14_insertBeg50,titre=u"Inserts distribution beginning of bacteria",NSI=NsiI_list(),ticks=5000,name="f12",save=True)

graph_multiline(AB_Pa14_insert[:,0:6000],ticks=200,line=6,name="f13",NSI=NsiI_list(),save=True)

AB_Pa14_insert_outNSI=group(AB_Pa14_insert[:,6011715:6111715],pas=10)
graph(AB_Pa14_insert_outNSI,titre=u"Inserts distribution of bacteria without NsiI site (by 10bp)",NSI=NsiI_list(),name="f13b",save=True)

print "# Table"

TablePics(AB_Pa14_insert,pic=2500)

YRS=motif()

fichier=open("fig/tableMotif.csv", "wb")
writer=csv.writer(fichier)
writer.writerow(("Position","Missmatch","A_max","A_mean","B_max","B_mean"))

for yrs in YRS:
    i,sens,missmatch=yrs
    table=[i,missmatch]
        
    if sens:
        for banque in [A_Pa14_insert, B_Pa14_insert]:
            table.append(np.max(banque[1:,i+49:i+81].sum(axis=0)))
            table.append(round(np.mean(banque[1:,i+49:i+81].sum(axis=0)),2))
            
    elif not sens:
        for banque in [A_Pa14_insert, B_Pa14_insert]:
            table.append(np.max(banque[1:,i-81:i-49].sum(axis=0)))
            table.append(round(np.mean(banque[1:,i-81:i-49].sum(axis=0)),2))
    writer.writerow(table)

fichier.close()

print "# End of encapsidation"

(Taille_A_Pa14,A_Pa14_insertEnd,A_Pa14_unused)=dataAnalysePaired(A_Pa14_DATA[0],A_Pa14_DATA[1],Size["Pa14"],Nsi=True)
histo(Taille_A_Pa14,titre="Histogram of inserts size",name="f14",axes=[0,750, 0,20000],save=True)
graph_global(A_Pa14_insertEnd,titre=u"End of A inserts distribution on PA14", name="f15",save=True)

print "# Left ends"

Ag_Pa14_DATA=csvData("DATA/A_Pa14_gauche.csv")
Ag_Pa14_insert=dataAnalyseStart(Ag_Pa14_DATA[0],Size["Pa14"])
(Taille_Ag_Pa14,Ag_Pa14_insertEnd,Ag_Pa14_unused)=dataAnalysePaired(Ag_Pa14_DATA[0],Ag_Pa14_DATA[1],Size["Pa14"],Nsi=True)

histo(Taille_Ag_Pa14,bins=range(8,max(Taille_Ag_Pa14),11),Xmax=350,name="f16",save=True)
graph_global(Ag_Pa14_insert,titre=u"A left inserts distribution on PA14", name="f17",save=True)

print "concentration closes to patterns = %f"%pourcentMotif(Ag_Pa14_insert,"Pa14")

print "# Without NsiI"

Ab30_Pa14_DATA=csvData("DATA/Ab30_Pa14_droite2.csv")
Ab30_Pa14_insert=dataAnalyseStart(Ab30_Pa14_DATA[0],Size["Pa14"])

print "concentration closes to patterns = %f"%pourcentMotif(Ab30_Pa14_insert,"Pa14")

print "# Others phages"

PM_DATA=csvData("DATA/PM105_PAOI_gauche.csv")
PM_insert=dataAnalyseStart(PM_DATA[0],Size["PAOI"])
PM=dataAnalysePaired(PM_DATA[0],PM_DATA[1],Size["PAOI"])
graph_global(PM_insert,titre=u"PM105 inserts distribution on PAOI", name="f20",save=True)
histo(PM[0],bins=range(8,max(PM[0]),11),Xmax=350,name="f22",save=True)

print "****"  
print "insertion moyenne =                %f"%np.mean(PM_insert[1:,:].sum(axis=0))
print "variance =                         %f"%np.var(PM_insert[1:,:].sum(axis=0))
print "concentration closes to patterns = %f"%pourcentMotif(PM_insert,"PAOI")

P1_DATA=csvData("DATA/2P1_II10_gauche.csv")
P1_insert=dataAnalyseStart(P1_DATA[0],Size["II10"])
P1=dataAnalysePaired(P1_DATA[0],P1_DATA[1],Size["II10"])
graph_global(P1_insert,titre=u"2P1 inserts distribution on II10", name="f18",save=True)
histo(P1[0],bins=range(8,max(P1[0]),11),Xmax=350,name="f19",save=True)

print "****"  
print "insertion moyenne =                %f"%np.mean(P1_insert[1:,:].sum(axis=0))
print "variance =                         %f"%np.var(P1_insert[1:,:].sum(axis=0))
print "concentration closes to patterns = %f"%pourcentMotif(P1_insert,"II10")

print "# ms223"

a_ms223=dataAnalyseStart(csvData("DATA/ms223_A.csv")[0],6500)
b_ms223=dataAnalyseStart(csvData("DATA/ms223_B.csv")[0],6500)
p1_ms223=dataAnalyseStart(csvData("DATA/ms223_2P1.csv")[0],6500)
pM_ms223=dataAnalyseStart(csvData("DATA/ms223_pm105.csv")[0],6500)
ab_ms223=addition([a_ms223,b_ms223])


A_ms223=np.array(a_ms223[1:,:].sum(axis=0)/a_ms223[1:,:].std(),dtype=float)
B_ms223=np.array(b_ms223[1:,:].sum(axis=0)/b_ms223[1:,:].std(),dtype=float)
AB_ms233=np.array(ab_ms223[1:,:].sum(axis=0)/ab_ms223[1:,:].std(),dtype=float)
P1_ms223=np.array(p1_ms223[1:,:].sum(axis=0)/p1_ms223[1:,:].std(),dtype=float)
PM_ms223=np.array(pM_ms223[1:,:].sum(axis=0)/pM_ms223[1:,:].std(),dtype=float)

plt.clf()
plt.plot(A_ms223)
plt.plot(B_ms223)
plt.plot(P1_ms223)
plt.plot(PM_ms223)

plt.legend((('A', 'B', '2P1', 'Pm105') ))
plt.xticks()
plt.savefig("fig/f23a.jpeg",dpi=600,transparant=True)


plt.clf()
plt.plot(A_ms223[3000:4000])
plt.plot(B_ms223[3000:4000])
plt.plot(P1_ms223[3000:4000])
plt.plot(PM_ms223[3000:4000])

plt.legend((('A', 'B', '2P1', 'Pm105') ))
plt.xticks()
plt.savefig("fig/f23b.jpeg",dpi=600,transparant=True)

t()

##March 2016

AB_Pa14_insert_outNSI=group(AB_Pa14_insert[:,6011715:6111715],pas=10) #Bigest  segment without NsiI
graph(AB_Pa14_insert_outNSI,titre=u"Inserts distribution of bacteria without NsiI site (by 10bp)",NSI=NsiI_list(),Ymax=150,name="c1",save=True) 
graphMirror(AB_Pa14_insert_outNSI,titre=u"Inserts distribution of bacteria without NsiI site (by 10bp)",NSI=NsiI_list(),Ymax=150,name="c2",save=True)

graphMirror(AB_Pa14_insert[:,6055600:6055700],NSI=NsiI_list(),name="c3",save=True) #degenerate hotspot

graphMirror(AB_Pa14_insert[:,5217624:5217700],titre=u"Inserts distribution at a hotspot position",ticks=25,name="c4",NSI=NsiI_list(),save=True) #f7-like

graphMirror(AB_Pa14_insert[:,5729500:5730000],titre=u"Inserts distribution at a hotspot position",ticks=100,name="c5",NSI=NsiI_list(),save=True) #poly hotspot

graphMirror(AB_Ab30_insert50,titre=u"A+B inserts distribution on Ab30 (group by 50bp)",name="c6",Ymax=500,NSI=range(6165,6170),save=True) #f4-like

AB_Pa14_insert1000=group(AB_Pa14_insert,1000)
graph(AB_Pa14_insert1000,titre=u"A+B inserts distribution on bacteria (group by 1000)",name="c7",Ymax=25000,save=True)
graphMirror(AB_Pa14_insert1000,titre=u"A+B inserts distribution on bacteria (group by 1000)",name="c8",Ymax=25000,save=True)

graphMirror(group(AB_Pa14_insert[:,7999:22000],10),name="c9",ticks=2000,Ymax=200,NSI=NsiI_list(),save=True)
graphMirror(group(AB_Pa14_insert[:,7999:22000],1),name="c10",ticks=2000,Ymax=50,NSI=NsiI_list(),save=True)

A_Pa14_insert254825=selectinrange(A_Pa14_DATA[0],range(254805,254830+1))

(Taille_A_Pa14_254825,A_Pa14_insertEnd_254825,A_Pa14_unused)=dataAnalysePaired(A_Pa14_insert254825,A_Pa14_DATA[1],Size["Pa14"],Nsi=True)
histoDouble(Taille_A_Pa14,Taille_A_Pa14_254825,titre="Histogram of inserts in hotspot ms211",name="n11",Xmax=750,save=True)

t()
raw_input("Done.")
