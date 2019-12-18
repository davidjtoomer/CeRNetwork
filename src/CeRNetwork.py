# import necessary python packages
import math
import pandas as pd

# CREATE THE FIRST PART OF CeRNETWORK: Mapping ceRNA-miRNA Interactions

# import the datasets for ceRNA:miRNA and delete extraneous columns
# TWO COLUMNS: 'ceRNA' and 'miRNA'
starBaseLnc = pd.read_csv('Databases/starBase_lncRNA.txt', delimiter = '\t')
starBaseSnc = pd.read_csv('Databases/starBase_sncRNA.txt', delimiter = '\t')
starBaseCirc = pd.read_csv('Databases/starBase_circRNA.txt', delimiter = '\t')
starBasePseudo = pd.read_csv('Databases/starBase_pseudogenes.txt', delimiter = '\t')
sbncIndexList = ['miRNAid', 'geneID', 'geneType', 'chromosome','start','end','strand','clipExpNum','degraExpNum','RBP','merClass','miRseq','align','targetSeq']
for index in sbncIndexList:
    del starBaseLnc[index]
    del starBaseSnc[index]
    del starBaseCirc[index]
    del starBasePseudo[index]
del starBaseLnc['pancancerNum']
del starBaseSnc['pancancerNum']
del starBasePseudo['pancancerNum']

# write the first part as a new csv file
ceRNetwork1 = starBaseLnc.append(starBaseSnc.append(starBaseCirc.append(starBasePseudo)), ignore_index = True) 
ceRNACRNList1 = []
miRNACRNList1 = []
for x in range(0, len(ceRNetwork1)):
    if ceRNetwork1['ceRNA'][x] not in ceRNACRNList1:
        ceRNACRNList1.append(ceRNetwork1['ceRNA'][x])
    if ceRNetwork1['miRNA'][x] not in miRNACRNList1:
        miRNACRNList1.append(ceRNetwork1['miRNA'][x])
print("CMN ceRNA Nodes:", len(ceRNACRNList1))
print("CMN miRNA Nodes:", len(miRNACRNList1))
ceRNetwork1Duplicates = ceRNetwork1.duplicated()
lenBeforeDrop1 = len(ceRNetwork1)
ceRNetwork1.drop_duplicates(subset = None, inplace = True)
print("CMN Edges:", len(ceRNetwork1))
ceRNetwork1.to_csv('ceRNetwork1.csv', index = False)

# CREATE THE SECOND PART OF CeRNETWORK: Mapping the miRNA-Gene Interactions

# import datasets for miRNA:mRNA and delete extraneous columns
# TWO COLUMNS: 'miRNA' and 'Gene'
miRTarBase = pd.read_csv('Databases/miRTarBase.tsv', delimiter = '\t')
mtbIndexList = ['miRTarBase ID', 'Species (miRNA)', 'Target Gene (Entrez Gene ID)','Species (Target Gene)','Experiments','Support Type','References (PMID)']
for index in mtbIndexList:
    del miRTarBase[index]

miRandaHighConserv = pd.read_csv('Databases/miRanda_HighScore_Conserved.txt', delimiter = '\t')
miRandaHighNonconserv = pd.read_csv('Databases/miRanda_HighScore_Nonconserved.txt', delimiter = '\t')
mrdIndexList = ['#mirbase_acc','gene_id','transcript_id','genBank_access','mirna_alignment','alignment','gene_alignment','mirna_start','mirna_end','gene_start','gene_end','genome_coordinates','conservation','align_score','seed_cat','energy','mirsvr_score']
for index in mrdIndexList:
    del miRandaHighConserv[index]
    del miRandaHighNonconserv[index]

miRecords = pd.read_csv('Databases/miRecords.tsv', delimiter = '\t')
mrcIndexList = ['pubmed_id','target_gene_species','genBank_access','miRNA_species']
for index in mrcIndexList:
    del miRecords[index]

tarBase = pd.read_csv('Databases/tarBase.csv', delimiter = '\t')
trbIndexList = ['geneId','species','cell_line','tissue','category','method','positive_negative','direct_indirect','up_down','condition']
nonHumanTargets = []
for x in range(0, len(tarBase)):
    if tarBase['species'][x] != "Homo sapiens":
        nonHumanTargets.append(x)
tarBase.drop(index = nonHumanTargets, inplace = True)
for index in trbIndexList:
    del tarBase[index]
    
# write the second part as a new csv file
ceRNetwork2 = miRTarBase.append(miRandaHighConserv.append(miRandaHighNonconserv.append(miRecords.append(tarBase))), ignore_index = True)
miRNACRNList2 = []
geneCRNList2 = []
for x in range(0, len(ceRNetwork2)):
    if ceRNetwork2['miRNA'][x] not in miRNACRNList2:
        miRNACRNList2.append(ceRNetwork2['miRNA'][x])
    if ceRNetwork2['Gene'][x] not in geneCRNList2:
        geneCRNList2.append(ceRNetwork2['Gene'][x])
print("MGN miRNA Nodes:", len(miRNACRNList2))
print("MGN Gene Nodes:", len(geneCRNList2))
ceRNetwork2Duplicates = ceRNetwork2.duplicated()
lenBeforeDrop2 = len(ceRNetwork2)
ceRNetwork2.drop_duplicates(subset = None, inplace = True)
print("MGN Edges:", len(ceRNetwork2))
ceRNetwork2.to_csv('ceRNetwork2.csv', index = False)

# CREATE THE THIRD PART OF CeRNETWORK: Mapping Gene-Disease Interactions

# import dataset for gene:disease networks
# TWO COLUMNS: 'Gene' and 'Disease'
disGeNet = pd.read_csv('Databases/disGeNet.tsv', delimiter = '\t')
dgnIndexList = ['geneId','diseaseName','score','NofPmids','NofSnps','source']
for index in dgnIndexList:
    del disGeNet[index]
    
# write the third part as a new csv file
ceRNetwork3 = disGeNet
ceRNetwork3.to_csv('ceRNetwork3.csv', index = False)

geneCRNList3 = []
diseaseCRNList3 = []
for x in range(0, len(ceRNetwork3)):
    if ceRNetwork3['Gene'][x] not in geneCRNList3:
        geneCRNList3.append(ceRNetwork3['Gene'][x])
    if ceRNetwork3['Disease'][x] not in diseaseCRNList3:
        diseaseCRNList3.append(ceRNetwork3['Disease'][x])
print("GDN Gene Nodes:", len(geneCRNList3))
print("GDN Disease Nodes:", len(diseaseCRNList3))
print("GDN Edges:", len(ceRNetwork3))

# CONSTRUCTING BIPARTITE NETWORKS: ceRNA-Gene Network (CGN)

print("Processing CGN1: Initializing miRNA Set")
cgnMiRList1 = []
for x in range(0, lenBeforeDrop2):
    if x == int(lenBeforeDrop2*(1/10)):
        print("10% complete...")
    if x == int(lenBeforeDrop2*(2/10)):
        print("20% complete...")
    if x == int(lenBeforeDrop2*(3/10)):
        print("30% complete...")
    if x == int(lenBeforeDrop2*(4/10)):
        print("40% complete...")
    if x == int(lenBeforeDrop2*(5/10)):
        print("50% complete...")
    if x == int(lenBeforeDrop2*(6/10)):
        print("60% complete...")
    if x == int(lenBeforeDrop2*(7/10)):
        print("70% complete...")
    if x == int(lenBeforeDrop2*(8/10)):
        print("80% complete...")
    if x == int(lenBeforeDrop2*(9/10)):
        print("90% complete...")
    if ceRNetwork2Duplicates[x] == False:
        cgnMiRList1.append(ceRNetwork2['miRNA'][x])
cgnMiRSet1 = set(cgnMiRList1)
print("miRNA Set 100% Initialized.")

print("Processing CGN1: Creating EdgeList 1")
ceRNACGN = []
miRNACGN = []
edgeListCGN1 = []
for x in range(0, lenBeforeDrop1):
    if x == int(lenBeforeDrop1*(1/10)):
        print("10% complete...")
    if x == int(lenBeforeDrop1*(2/10)):
        print("20% complete...")
    if x == int(lenBeforeDrop1*(3/10)):
        print("30% complete...")
    if x == int(lenBeforeDrop1*(4/10)):
        print("40% complete...")
    if x == int(lenBeforeDrop1*(5/10)):
        print("50% complete...")
    if x == int(lenBeforeDrop1*(6/10)):
        print("60% complete...")
    if x == int(lenBeforeDrop1*(7/10)):
        print("70% complete...")
    if x == int(lenBeforeDrop1*(8/10)):
        print("80% complete...")
    if x == int(lenBeforeDrop1*(9/10)):
        print("90% complete...")
    if ceRNetwork1Duplicates[x] == False:
        if ceRNetwork1['miRNA'][x] in cgnMiRSet1:
            ceRNACGN.append(ceRNetwork1['ceRNA'][x])
            miRNACGN.append(ceRNetwork1['miRNA'][x])
            edgeListCGN1.append([ceRNetwork1['ceRNA'][x], ceRNetwork1['miRNA'][x]])
print("CGN1 EdgeList 100% Complete.")

print("Processing CGN2: Initializing miRNA Set")
cgnMiRList2 = []
for x in range(0, lenBeforeDrop1):
    if ceRNetwork1Duplicates[x] == False:
        cgnMiRList2.append(ceRNetwork1['miRNA'][x])
cgnMiRSet2 = set(cgnMiRList2)
print("miRNA Set 100% Initialized.")

print("Processing CGN2: Creating EdgeList 2")
geneCGN = []
edgeListCGN2 = []
for x in range(0, lenBeforeDrop2):
    if x == int(lenBeforeDrop2*(1/10)):
        print("10% complete...")
    if x == int(lenBeforeDrop2*(2/10)):
        print("20% complete...")
    if x == int(lenBeforeDrop2*(3/10)):
        print("30% complete...")
    if x == int(lenBeforeDrop2*(4/10)):
        print("40% complete...")
    if x == int(lenBeforeDrop2*(5/10)):
        print("50% complete...")
    if x == int(lenBeforeDrop2*(6/10)):
        print("60% complete...")
    if x == int(lenBeforeDrop2*(7/10)):
        print("70% complete...")
    if x == int(lenBeforeDrop2*(8/10)):
        print("80% complete...")
    if x == int(lenBeforeDrop2*(9/10)):
        print("90% complete...")
    if ceRNetwork2Duplicates[x] == False:
        if ceRNetwork2['miRNA'][x] in cgnMiRSet2:
            geneCGN.append(ceRNetwork2['Gene'][x])
            edgeListCGN2.append([ceRNetwork2['miRNA'][x], ceRNetwork2['Gene'][x]])
print("CGN2 EdgeList 100% Complete.")

# method for dropping duplicate items in a list
def drop_duplicates(List):
    set = []
    for line in List:
        if line not in set:
            set.append(line)
    return set

# find number of nodes and edges
ceRListCGN = drop_duplicates(ceRNACGN)
lenCeRListCGN = len(ceRListCGN)
print("CMGN ceRNA Nodes:", lenCeRListCGN)

miRListCGN = drop_duplicates(miRNACGN)
lenMiRListCGN = len(miRListCGN)
print("CMGN miRNA Nodes:", lenMiRListCGN)

geneListCGN = drop_duplicates(geneCGN)
lenGeneListCGN = len(geneListCGN)
print("CMGN Gene Nodes:", lenGeneListCGN)

numEdgesCGN1 = len(edgeListCGN1)
numEdgesCGN2 = len(edgeListCGN2)
numEdgesCGN = numEdgesCGN1 + numEdgesCGN2
print("CMGN Edges:", numEdgesCGN)

print("Initializing Hypergeometric Test...")
# find average degrees of the bipartite networks
avgDegCGN1 = (2*numEdgesCGN1)/(lenCeRListCGN + lenMiRListCGN)
avgDegCGN2 = (2*numEdgesCGN2)/(lenMiRListCGN + lenGeneListCGN)

# hypergeometric test: statistical features of a random network 
meanCGN = (avgDegCGN1*avgDegCGN2)/lenMiRListCGN
sigmaCGN = math.sqrt(meanCGN*(1-(avgDegCGN1/lenMiRListCGN)))
print("CGN Mean:", meanCGN)
print("CGN Standard Deviation:", sigmaCGN)

# returns the union of two sets
def union(List1, List2):
    set = []
    for x in range(0, len(List1)):
        if List1[x] in List2:
            set.append(List1[x])
    return set

# create master lists of miRNAs associated with each ceRNA or gene node
superMiRListCGN1 = []
print("Processing: Creating the first large boy list of miRNAs")
for x in range(len(ceRListCGN)):
    if x == int(len(ceRListCGN)*(1/10)):
        print("10% complete...")
    if x == int(len(ceRListCGN)*(2/10)):
        print("20% complete...")
    if x == int(len(ceRListCGN)*(3/10)):
        print("30% complete...")
    if x == int(len(ceRListCGN)*(4/10)):
        print("40% complete...")
    if x == int(len(ceRListCGN)*(5/10)):
        print("50% complete...")
    if x == int(len(ceRListCGN)*(6/10)):
        print("60% complete...")
    if x == int(len(ceRListCGN)*(7/10)):
        print("70% complete...")
    if x == int(len(ceRListCGN)*(8/10)):
        print("80% complete...")
    if x == int(len(ceRListCGN)*(9/10)):
        print("90% complete...")
    mirsetCGN1 = [edgeListCGN1[y][1] for y in range(len(edgeListCGN1)) if ceRListCGN[x] == edgeListCGN1[y][0]]
    superMiRListCGN1.append(mirsetCGN1)
print("First Master miRNA List 100% Complete.")

print("Processing: Creating the second large boy list of miRNAs")
superMiRListCGN2 = []
for x in range(len(geneListCGN)):
    if x == int(len(geneListCGN)*(1/10)):
        print("10% complete...")
    if x == int(len(geneListCGN)*(2/10)):
        print("20% complete...")
    if x == int(len(geneListCGN)*(3/10)):
        print("30% complete...")
    if x == int(len(geneListCGN)*(4/10)):
        print("40% complete...")
    if x == int(len(geneListCGN)*(5/10)):
        print("50% complete...")
    if x == int(len(geneListCGN)*(6/10)):
        print("60% complete...")
    if x == int(len(geneListCGN)*(7/10)):
        print("70% complete...")
    if x == int(len(geneListCGN)*(8/10)):
        print("80% complete...")
    if x == int(len(geneListCGN)*(9/10)):
        print("90% complete...")
    mirsetCGN2 = [edgeListCGN2[y][0] for y in range(len(edgeListCGN2)) if geneListCGN == edgeListCGN2[y][1]]
    superMiRListCGN2.append(mirsetCGN2)
print("Second Master miRNA List 100% Complete.")

# find the unions of each subset of the master list
print("Finding unions of the miRNA lists:")
CGN = []
for x in range(len(ceRListCGN)):
    if x == int(len(ceRListCGN)*(1/10)):
        print("10% complete...")
    if x == int(len(ceRListCGN)*(2/10)):
        print("20% complete...")
    if x == int(len(ceRListCGN)*(3/10)):
        print("30% complete...")
    if x == int(len(ceRListCGN)*(4/10)):
        print("40% complete...")
    if x == int(len(ceRListCGN)*(5/10)):
        print("50% complete...")
    if x == int(len(ceRListCGN)*(6/10)):
        print("60% complete...")
    if x == int(len(ceRListCGN)*(7/10)):
        print("70% complete...")
    if x == int(len(ceRListCGN)*(8/10)):
        print("80% complete...")
    if x == int(len(ceRListCGN)*(9/10)):
        print("90% complete...")
    list1 = superMiRListCGN1[x]
    for y in range(len(geneListCGN)):
        list2 = superMiRListCGN2[y]
        count = len(union(list1, list2))
        zScore = (count - meanCGN)/sigmaCGN
        if zScore > 3:
            CGN.append([ceRListCGN[x], geneListCGN[y], zScore])
print("CGN 100% Complete.")
print("Creating node and edge lists...")
# create node and edge lists for the CGN after hypergeometric reduction
ceRNACGNList = []
geneCGNList = []
for x in range(0, len(CGN)):
    if CGN[x][0] not in ceRNACGNList:
        ceRNACGNList.append(CGN[x][0])
    if CGN[x][1] not in geneCGNList:
        geneCGNList.append(CGN[x][1])
print("CGN ceRNA Nodes:", len(ceRNACGNList))
print("CGN Gene Nodes:", len(geneCGNList))
print("CGN Edges:", len(CGN))

# write the list as a pandas dataframe and export as a csv file
columnsCGN = ['ceRNA','Gene','Score']
dfCGN = pd.DataFrame(CGN, columns = columnsCGN)
dfCGN.to_csv('ceRNetwork_CMN.csv', index = False)

print("Processing: Creating the First Master miRNA List")
superMiRListCGN1 = []
for x in range(len(ceRListCGN)):
    if x == int(len(ceRListCGN)*(1/10)):
        print("10% complete...")
    if x == int(len(ceRListCGN)*(2/10)):
        print("20% complete...")
    if x == int(len(ceRListCGN)*(3/10)):
        print("30% complete...")
    if x == int(len(ceRListCGN)*(4/10)):
        print("40% complete...")
    if x == int(len(ceRListCGN)*(5/10)):
        print("50% complete...")
    if x == int(len(ceRListCGN)*(6/10)):
        print("60% complete...")
    if x == int(len(ceRListCGN)*(7/10)):
        print("70% complete...")
    if x == int(len(ceRListCGN)*(8/10)):
        print("80% complete...")
    if x == int(len(ceRListCGN)*(9/10)):
        print("90% complete...")
    mirsetCGN1 = []
    for y in range(len(edgeListCGN1)):
        if ceRListCGN[x] == edgeListCGN1[y][0]:
            mirsetCGN1.append(edgeListCGN1[y][1])
    superMiRListCGN1.append(mirsetCGN1)
print("First miRNA list complete.")

print("Processing: Creating the Second Master miRNA List")
superMiRListCGN2 = []
for x in range(len(geneListCGN)):
    if x == int(len(geneListCGN)*(1/10)):
        print("10% complete...")
    if x == int(len(geneListCGN)*(2/10)):
        print("20% complete...")
    if x == int(len(geneListCGN)*(3/10)):
        print("30% complete...")
    if x == int(len(geneListCGN)*(4/10)):
        print("40% complete...")
    if x == int(len(geneListCGN)*(5/10)):
        print("50% complete...")
    if x == int(len(geneListCGN)*(6/10)):
        print("60% complete...")
    if x == int(len(geneListCGN)*(7/10)):
        print("70% complete...")
    if x == int(len(geneListCGN)*(8/10)):
        print("80% complete...")
    if x == int(len(geneListCGN)*(9/10)):
        print("90% complete...")
    mirsetCGN2 = []
    for y in range(len(edgeListCGN2)):
        if geneListCGN[x] == edgeListCGN2[y][1]:
            mirsetCGN2.append(edgeListCGN2[y][0])
    superMiRListCGN2.append(mirsetCGN2)
print("Second miRNA list complete.")

print("Finding unions of the miRNA lists:")
CGN = []
for x in range(len(ceRListCGN)):
    if x == int(len(ceRListCGN)*(1/10)):
        print("10% complete...")
    if x == int(len(ceRListCGN)*(2/10)):
        print("20% complete...")
    if x == int(len(ceRListCGN)*(3/10)):
        print("30% complete...")
    if x == int(len(ceRListCGN)*(4/10)):
        print("40% complete...")
    if x == int(len(ceRListCGN)*(5/10)):
        print("50% complete...")
    if x == int(len(ceRListCGN)*(6/10)):
        print("60% complete...")
    if x == int(len(ceRListCGN)*(7/10)):
        print("70% complete...")
    if x == int(len(ceRListCGN)*(8/10)):
        print("80% complete...")
    if x == int(len(ceRListCGN)*(9/10)):
        print("90% complete...")
    list1 = superMiRListCGN1[x]
    for y in range(len(geneListCGN)):
        list2 = superMiRListCGN2[y]
        count = len(union(list1, list2))
        zScore = (count - meanCGN)/sigmaCGN
        if zScore > 3:
            CGN.append([ceRListCGN[x], geneListCGN[y], zScore])
print("CGN 100% Complete.")
print("Creating node and edge lists...")
# create node and edge lists for the CGN after hypergeometric reduction
ceRNACGNList = []
geneCGNList = []
for x in range(0, len(CGN)):
    if CGN[x][0] not in ceRNACGNList:
        ceRNACGNList.append(CGN[x][0])
    if CGN[x][1] not in geneCGNList:
        geneCGNList.append(CGN[x][1])
print("CGN ceRNA Nodes:", len(ceRNACGNList))
print("CGN Gene Nodes:", len(geneCGNList))
print("CGN Edges:", len(CGN))

# write the list as a pandas dataframe and export as a csv file
columnsCGN = ['ceRNA','Gene','Score']
dfCGN = pd.DataFrame(CGN, columns = columnsCGN)
dfCGN.to_csv('ceRNetwork_CMN.csv', index = False)

# CONSTRUCTING BIPARTITE NETWORKS: miRNA-Disease Network (MDN)

# find the number of miRNAs, genes, and diseases in the network (combination of CRN2 and CRN3)
print("Processing MDN1: Initializing Gene Set 1")
mdnGeneList1 = []
for x in range(0, len(ceRNetwork3)):
    if x == int(len(ceRNetwork3)*(1/10)):
        print("10% complete...")
    if x == int(len(ceRNetwork3)*(2/10)):
        print("20% complete...")
    if x == int(len(ceRNetwork3)*(3/10)):
        print("30% complete...")
    if x == int(len(ceRNetwork3)*(4/10)):
        print("40% complete...")
    if x == int(len(ceRNetwork3)*(5/10)):
        print("50% complete...")
    if x == int(len(ceRNetwork3)*(6/10)):
        print("60% complete...")
    if x == int(len(ceRNetwork3)*(7/10)):
        print("70% complete...")
    if x == int(len(ceRNetwork3)*(8/10)):
        print("80% complete...")
    if x == int(len(ceRNetwork3)*(9/10)):
        print("90% complete...")
    mdnGeneList1.append(ceRNetwork3['Gene'][x])
mdnGeneSet1 = set(mdnGeneList1)
print("Gene Set 1 100% Initialized.")

print("Processing MDN1: Creating EdgeList 1")
miRNAMDN = []
geneMDN = []
edgeListMDN1 = []
for x in range(0, lenBeforeDrop2):
    if x == int(lenBeforeDrop2*(1/10)):
        print("10% complete...")
    if x == int(lenBeforeDrop2*(2/10)):
        print("20% complete...")
    if x == int(lenBeforeDrop2*(3/10)):
        print("30% complete...")
    if x == int(lenBeforeDrop2*(4/10)):
        print("40% complete...")
    if x == int(lenBeforeDrop2*(5/10)):
        print("50% complete...")
    if x == int(lenBeforeDrop2*(6/10)):
        print("60% complete...")
    if x == int(lenBeforeDrop2*(7/10)):
        print("70% complete...")
    if x == int(lenBeforeDrop2*(8/10)):
        print("80% complete...")
    if x == int(lenBeforeDrop2*(9/10)):
        print("90% complete...")
    if ceRNetwork2Duplicates[x] == False:
        if ceRNetwork2['Gene'][x] in mdnGeneSet1:
            miRNAMDN.append(ceRNetwork2['miRNA'][x])
            geneMDN.append(ceRNetwork2['Gene'][x])
            edgeListMDN1.append([ceRNetwork2['miRNA'][x], ceRNetwork2['Gene'][x]])
print("MDN1 EdgeList 100% Complete.")

print("Processing MDN2: Initializing Gene Set 2")
mdnGeneList2 = []
for x in range(0, lenBeforeDrop2):
    if x == int(lenBeforeDrop2*(1/10)):
        print("10% complete...")
    if x == int(lenBeforeDrop2*(2/10)):
        print("20% complete...")
    if x == int(lenBeforeDrop2*(3/10)):
        print("30% complete...")
    if x == int(lenBeforeDrop2*(4/10)):
        print("40% complete...")
    if x == int(lenBeforeDrop2*(5/10)):
        print("50% complete...")
    if x == int(lenBeforeDrop2*(6/10)):
        print("60% complete...")
    if x == int(lenBeforeDrop2*(7/10)):
        print("70% complete...")
    if x == int(lenBeforeDrop2*(8/10)):
        print("80% complete...")
    if x == int(lenBeforeDrop2*(9/10)):
        print("90% complete...")
    if ceRNetwork2Duplicates[x] == False:
        mdnGeneList2.append(ceRNetwork2['Gene'][x])
mdnGeneSet2 = set(mdnGeneList2)
print("Gene Set 2 100% Initialized.")

print("Processing MDN2: Creating EdgeList 2")
diseaseMDN = []
edgeListMDN2 = []
for x in range(0, len(ceRNetwork3)):
    if x == int(len(ceRNetwork3)*(1/10)):
        print("10% complete...")
    if x == int(len(ceRNetwork3)*(2/10)):
        print("20% complete...")
    if x == int(len(ceRNetwork3)*(3/10)):
        print("30% complete...")
    if x == int(len(ceRNetwork3)*(4/10)):
        print("40% complete...")
    if x == int(len(ceRNetwork3)*(5/10)):
        print("50% complete...")
    if x == int(len(ceRNetwork3)*(6/10)):
        print("60% complete...")
    if x == int(len(ceRNetwork3)*(7/10)):
        print("70% complete...")
    if x == int(len(ceRNetwork3)*(8/10)):
        print("80% complete...")
    if x == int(len(ceRNetwork3)*(9/10)):
        print("90% complete...")
    if ceRNetwork3['Gene'][x] in mdnGeneSet2:
        diseaseMDN.append(ceRNetwork3['Disease'][x])
        edgeListMDN2.append([ceRNetwork3['Gene'][x], ceRNetwork3['Disease'][x]])
print("MDN2 100% Complete.")

# find number of nodes and edges
miRListMDN = drop_duplicates(miRNAMDN)
lenMiRListMDN = len(miRListMDN)
print("MGDN miRNA Nodes:", lenMiRListMDN)

geneListMDN = drop_duplicates(geneMDN)
lenGeneListMDN1 = len(geneListMDN)

print("MGDN Gene Nodes:", lenGeneListMDN1)

diseaseListMDN = drop_duplicates(diseaseMDN)
lenDiseaseListMDN = len(diseaseListMDN)
print("MGDN Disease Nodes:", lenDiseaseListMDN)

numEdgesMDN1 = len(edgeListMDN1)
numEdgesMDN2 = len(edgeListMDN2)
numEdgesMDN = numEdgesMDN1 + numEdgesMDN2
print("MGDN Edges:", numEdgesMDN)

print("Initializing Hypergeometric Test")
# find average degrees of the bipartite networks
avgDegMDN1 = (2*numEdgesMDN1)/(lenMiRListMDN + lenGeneListMDN1)
avgDegMDN2 = (2*numEdgesMDN2)/(lenGeneListMDN1 + lenDiseaseListMDN)

# hypergeometric test
meanMDN = (avgDegMDN1*avgDegMDN2)/lenGeneListMDN1
sigmaMDN = math.sqrt(meanMDN*(1-(avgDegMDN1/lenGeneListMDN1)))
print("MDN Mean:", meanMDN)
print("MDN Standard Deviation:", sigmaMDN)

# create universal gene lists for each miRNA and disease node
superGeneListMDN1 = []
print("Processing: Creating the First Master Gene List")
for x in range(len(miRListMDN)):
    if x == int(len(miRListMDN)*(1/10)):
        print("10% complete...")
    if x == int(len(miRListMDN)*(2/10)):
        print("20% complete...")
    if x == int(len(miRListMDN)*(3/10)):
        print("30% complete...")
    if x == int(len(miRListMDN)*(4/10)):
        print("40% complete...")
    if x == int(len(miRListMDN)*(5/10)):
        print("50% complete...")
    if x == int(len(miRListMDN)*(6/10)):
        print("60% complete...")
    if x == int(len(miRListMDN)*(7/10)):
        print("70% complete...")
    if x == int(len(miRListMDN)*(8/10)):
        print("80% complete...")
    if x == int(len(miRListMDN)*(9/10)):
        print("90% complete...")
    genesetMDN1 = []
    for y in range(len(edgeListMDN1)):
        if miRListMDN[x] == edgeListMDN1[y][0]:
            genesetMDN1.append(edgeListMDN1[y][1])
    superGeneListMDN1.append(genesetMDN1)
print("First Master miRNA List 100% Complete.")

print("MGDN Gene Nodes:", lenGeneListMDN1)

diseaseListMDN = drop_duplicates(diseaseMDN)
lenDiseaseListMDN = len(diseaseListMDN)
print("MGDN Disease Nodes:", lenDiseaseListMDN)

numEdgesMDN1 = len(edgeListMDN1)
numEdgesMDN2 = len(edgeListMDN2)
numEdgesMDN = numEdgesMDN1 + numEdgesMDN2
print("MGDN Edges:", numEdgesMDN)

print("Initializing Hypergeometric Test")
# find average degrees of the bipartite networks
avgDegMDN1 = (2*numEdgesMDN1)/(lenMiRListMDN + lenGeneListMDN1)
avgDegMDN2 = (2*numEdgesMDN2)/(lenGeneListMDN1 + lenDiseaseListMDN)

# hypergeometric test
meanMDN = (avgDegMDN1*avgDegMDN2)/lenGeneListMDN1
sigmaMDN = math.sqrt(meanMDN*(1-(avgDegMDN1/lenGeneListMDN1)))
print("MDN Mean:", meanMDN)
print("MDN Standard Deviation:", sigmaMDN)

# create universal gene lists for each miRNA and disease node
superGeneListMDN1 = []
print("Processing: Creating the First Master Gene List")
for x in range(len(miRListMDN)):
    if x == int(len(miRListMDN)*(1/10)):
        print("10% complete...")
    if x == int(len(miRListMDN)*(2/10)):
        print("20% complete...")
    if x == int(len(miRListMDN)*(3/10)):
        print("30% complete...")
    if x == int(len(miRListMDN)*(4/10)):
        print("40% complete...")
    if x == int(len(miRListMDN)*(5/10)):
        print("50% complete...")
    if x == int(len(miRListMDN)*(6/10)):
        print("60% complete...")
    if x == int(len(miRListMDN)*(7/10)):
        print("70% complete...")
    if x == int(len(miRListMDN)*(8/10)):
        print("80% complete...")
    if x == int(len(miRListMDN)*(9/10)):
        print("90% complete...")
    genesetMDN1 = []
    for y in range(len(edgeListMDN1)):
        if miRListMDN[x] == edgeListMDN1[y][0]:
            genesetMDN1.append(edgeListMDN1[y][1])
    superGeneListMDN1.append(genesetMDN1)
print("First Master miRNA List 100% Complete.")

print("Finding Unions of the Master Gene Lists.")
MDN = []
for x in range(len(miRListMDN)):
    if x == int(len(miRListMDN)*(1/10)):
        print("10% complete...")
    if x == int(len(miRListMDN)*(2/10)):
        print("20% complete...")
    if x == int(len(miRListMDN)*(3/10)):
        print("30% complete...")
    if x == int(len(miRListMDN)*(4/10)):
        print("40% complete...")
    if x == int(len(miRListMDN)*(5/10)):
        print("50% complete...")
    if x == int(len(miRListMDN)*(6/10)):
        print("60% complete...")
    if x == int(len(miRListMDN)*(7/10)):
        print("70% complete...")
    if x == int(len(miRListMDN)*(8/10)):
        print("80% complete...")
    if x == int(len(miRListMDN)*(9/10)):
        print("90% complete...")
    list1 = superGeneListMDN1[x]
    for y in range(len(diseaseListMDN)):
        list2 = superGeneListMDN2[y]
        count = len(union(list1, list2))
        zScore = (count - meanMDN)/sigmaMDN
        if zScore > 3:
            MDN.append([miRListMDN[x], diseaseListMDN[y], zScore])
print("MDN 100% Complete.")

# create node and edge lists for the CGN after hypergeometric reduction
print("Creating node and edge lists...")
miRNAMDNList = []
diseaseMDNList = []
for x in range(0, len(MDN)):
    if MDN[x][0] not in miRNAMDNList:
        miRNAMDNList.append(MDN[x][0])
    if MDN[x][1] not in diseaseMDNList:
        diseaseMDNList.append(MDN[x][1])
print("MDN miRNA Nodes:", len(miRNAMDNList))
print("MDN Disease Nodes:", len(diseaseMDNList))
print("MDN Edges:", len(MDN))

# write the list as a pandas dataframe and export as a csv file
columnsMDN = ['miRNA', 'Disease','Score']
dfMDN = pd.DataFrame(MDN, columns = columnsMDN)
dfMDN.to_csv('ceRNetwork_MDN.csv', index = False)

# CONSTRUCTING BIPARTITE NETWORKS: ceRNA-Disease Network (CDN)

# find the number of ceRNAs, miRNAs, and diseases in the network (combination of CRN1 and MDN)
print("Processing CDN1: Initializing miRNA Set 1")
cdnMiRList1 = []
for x in range(len(dfMDN)):
    if x == int(len(dfMDN)*(1/10)):
        print("10% complete...")
    if x == int(len(dfMDN)*(2/10)):
        print("20% complete...")
    if x == int(len(dfMDN)*(3/10)):
        print("30% complete...")
    if x == int(len(dfMDN)*(4/10)):
        print("40% complete...")
    if x == int(len(dfMDN)*(5/10)):
        print("50% complete...")
    if x == int(len(dfMDN)*(6/10)):
        print("60% complete...")
    if x == int(len(dfMDN)*(7/10)):
        print("70% complete...")
    if x == int(len(dfMDN)*(8/10)):
        print("80% complete...")
    if x == int(len(dfMDN)*(9/10)):
        print("90% complete...")
    cdnMiRList1.append(dfMDN['miRNA'][x])
cdnMiRSet1 = set(cdnMiRList1)
print("miRNA Set 1 100% Initialized.")

print("Processing CDN1: Creating EdgeList 1")
ceRNACDN = []
miRNACDN = []
edgeListCDN1 = []
for x in range(0, lenBeforeDrop1):
    if x == int(lenBeforeDrop1*(1/10)):
        print("10% complete...")
    if x == int(lenBeforeDrop1*(2/10)):
        print("20% complete...")
    if x == int(lenBeforeDrop1*(3/10)):
        print("30% complete...")
    if x == int(lenBeforeDrop1*(4/10)):
        print("40% complete...")
    if x == int(lenBeforeDrop1*(5/10)):
        print("50% complete...")
    if x == int(lenBeforeDrop1*(6/10)):
        print("60% complete...")
    if x == int(lenBeforeDrop1*(7/10)):
        print("70% complete...")
    if x == int(lenBeforeDrop1*(8/10)):
        print("80% complete...")
    if x == int(lenBeforeDrop1*(9/10)):
        print("90% complete...")
    if ceRNetwork1Duplicates[x] == False:
        if ceRNetwork1['miRNA'][x] in cdnMiRSet1:
            ceRNACDN.append(ceRNetwork1['ceRNA'][x])
            miRNACDN.append(ceRNetwork1['miRNA'][x])
            edgeListCDN1.append([ceRNetwork1['ceRNA'][x], ceRNetwork1['miRNA'][x]])
print("EdgeList 1 100% Complete.")

print("Processing CDN2: Initializing miRNA Set 2")
cdnMiRList2 = []
cdnMiRSet2 = cgnMiRSet2
print("miRNA Set 2 100% Initialized.")

print("Processing CDN2: Creating EdgeList 2")
diseaseCDN = []
edgeListCDN2 = []
for x in range(0, len(dfMDN)):
    if x == int(len(dfMDN)*(1/10)):
        print("10% complete...")
    if x == int(len(dfMDN)*(2/10)):
        print("20% complete...")
    if x == int(len(dfMDN)*(3/10)):
        print("30% complete...")
    if x == int(len(dfMDN)*(4/10)):
        print("40% complete...")
    if x == int(len(dfMDN)*(5/10)):
        print("50% complete...")
    if x == int(len(dfMDN)*(6/10)):
        print("60% complete...")
    if x == int(len(dfMDN)*(7/10)):
        print("70% complete...")
    if x == int(len(dfMDN)*(8/10)):
        print("80% complete...")
    if x == int(len(dfMDN)*(9/10)):
        print("90% complete...")
    if dfMDN['miRNA'][x] in cdnMiRSet2:
        diseaseCDN.append(dfMDN['Disease'][x])
        edgeListCDN2.append([dfMDN['miRNA'][x], dfMDN['Disease'][x]])
print("EdgeList 2 100% Complete.")

# find number of nodes and edges
ceRListCDN = drop_duplicates(ceRNACDN)
lenCeRListCDN = len(ceRListCDN)
print("CMDN ceRNA Nodes:", lenCeRListCDN)

miRListCDN = drop_duplicates(miRNACDN)
lenMiRListCDN = len(miRListCDN)
print("CMDN miRNA Nodes:", lenMiRListCDN)

diseaseListCDN = drop_duplicates(diseaseCDN)
lenDiseaseListCDN = len(diseaseListCDN)
print("CMDN Disease Nodes:", lenDiseaseListCDN)

numEdgesCDN1 = len(edgeListCDN1)
numEdgesCDN2 = len(edgeListCDN2)
numEdgesCDN = numEdgesCDN1 + numEdgesCDN2
print("CMDN Edges:", numEdgesCDN)

print("Initializing Hypergeometric Test...")
# find average degrees of the bipartite networks
avgDegCDN1 = (2*numEdgesCDN1)/(lenCeRListCDN + lenMiRListCDN)
avgDegCDN2 = (2*numEdgesCDN2)/(lenMiRListCDN + lenDiseaseListCDN)

# hypergeometric test
meanCDN = (avgDegCDN1*avgDegCDN2)/lenMiRListCDN
sigmaCDN = math.sqrt(meanCDN*(1-(avgDegCDN1/lenMiRListCDN)))
print("CDN Mean:", meanCDN)
print("CDN Standard Deviation:", sigmaCDN)

print("Processing: Creating the First Master miRNA List")
superMiRListCDN1 = []
for x in range(len(ceRListCDN)):
    if x == int(len(ceRListCDN)*(1/10)):
        print("10% complete...")
    if x == int(len(ceRListCDN)*(2/10)):
        print("20% complete...")
    if x == int(len(ceRListCDN)*(3/10)):
        print("30% complete...")
    if x == int(len(ceRListCDN)*(4/10)):
        print("40% complete...")
    if x == int(len(ceRListCDN)*(5/10)):
        print("50% complete...")
    if x == int(len(ceRListCDN)*(6/10)):
        print("60% complete...")
    if x == int(len(ceRListCDN)*(7/10)):
        print("70% complete...")
    if x == int(len(ceRListCDN)*(8/10)):
        print("80% complete...")
    if x == int(len(ceRListCDN)*(9/10)):
        print("90% complete...")
    mirsetCDN1 = []
    for y in range(len(edgeListCDN1)):
        if ceRListCDN[x] == edgeListCDN1[y][0]:
            mirsetCDN1.append(edgeListCDN1[y][1])
    superMiRListCDN1.append(mirsetCDN1)
print("First Master miRNA list complete.")

print("Processing: Creating the Second Master miRNA List")
superMiRListCDN2 = []
for x in range(len(diseaseListCDN)):
    if x == int(len(diseaseListCDN)*(1/10)):
        print("10% complete...")
    if x == int(len(diseaseListCDN)*(2/10)):
        print("20% complete...")
    if x == int(len(diseaseListCDN)*(3/10)):
        print("30% complete...")
    if x == int(len(diseaseListCDN)*(4/10)):
        print("40% complete...")
    if x == int(len(diseaseListCDN)*(5/10)):
        print("50% complete...")
    if x == int(len(diseaseListCDN)*(6/10)):
        print("60% complete...")
    if x == int(len(diseaseListCDN)*(7/10)):
        print("70% complete...")
    if x == int(len(diseaseListCDN)*(8/10)):
        print("80% complete...")
    if x == int(len(diseaseListCDN)*(9/10)):
        print("90% complete...")
    mirsetCDN2 = []
    for y in range(len(edgeListCDN2)):
        if diseaseListCDN[x] == edgeListCDN2[y][1]:
            mirsetCDN2.append(edgeListCDN2[y][0])
    superMiRListCDN2.append(mirsetCDN2)
print("Second Master miRNA list complete.")

print("Creating CDN: Finding intersections of the miRNA lists")
CDN = []
for x in range(len(ceRListCDN)):
    if x == int(len(ceRListCDN)*(1/10)):
        print("10% complete...")
    if x == int(len(ceRListCDN)*(2/10)):
        print("20% complete...")
    if x == int(len(ceRListCDN)*(3/10)):
        print("30% complete...")
    if x == int(len(ceRListCDN)*(4/10)):
        print("40% complete...")
    if x == int(len(ceRListCDN)*(5/10)):
        print("50% complete...")
    if x == int(len(ceRListCDN)*(6/10)):
        print("60% complete...")
    if x == int(len(ceRListCDN)*(7/10)):
        print("70% complete...")
    if x == int(len(ceRListCDN)*(8/10)):
        print("80% complete...")
    if x == int(len(ceRListCDN)*(9/10)):
        print("90% complete...")
    list1 = superMiRListCDN1[x]
    for y in range(len(diseaseListCDN)):
        list2 = superMiRListCDN2[y]
        count = len(union(list1, list2))
        zScore = (count - meanCDN)/sigmaCDN
        if zScore > 3:
            CDN.append([ceRListCDN[x], diseaseListCDN[y], zScore])
print("CDN 100% Complete.")

# create node and edge lists for the CGN after hypergeometric reduction
print("Creating node and edge lists...")
ceRNACDNList = []
diseaseCDNList = []
for x in range(0, len(CDN)):
    if CDN[x][0] not in ceRNACDNList:
        ceRNACDNList.append(CDN[x][0])
    if CDN[x][1] not in diseaseCDNList:
        diseaseCDNList.append(CDN[x][1])
print("CDN ceRNA Nodes:", len(ceRNACDNList))
print("CDN Disease Nodes:", len(diseaseCDNList))
print("CDN Edges:", len(CDN))

# write the list as a pandas dataframe and export as a csv file
columnsCDN = ['ceRNA','Disease','Score']
dfCGN = pd.DataFrame(CDN, columns = columnsCDN)
dfCGN.to_csv('ceRNetwork_CDN.csv', index = False)

print("Processing: Creating the Second Master miRNA List")
superMiRListCDN2 = []
for x in range(len(diseaseListCDN)):
    if x == int(len(diseaseListCDN)*(1/10)):
        print("10% complete...")
    if x == int(len(diseaseListCDN)*(2/10)):
        print("20% complete...")
    if x == int(len(diseaseListCDN)*(3/10)):
        print("30% complete...")
    if x == int(len(diseaseListCDN)*(4/10)):
        print("40% complete...")
    if x == int(len(diseaseListCDN)*(5/10)):
        print("50% complete...")
    if x == int(len(diseaseListCDN)*(6/10)):
        print("60% complete...")
    if x == int(len(diseaseListCDN)*(7/10)):
        print("70% complete...")
    if x == int(len(diseaseListCDN)*(8/10)):
        print("80% complete...")
    if x == int(len(diseaseListCDN)*(9/10)):
        print("90% complete...")
    mirsetCDN2 = []
    for y in range(len(edgeListCDN2)):
        if diseaseListCDN[x] == edgeListCDN2[y][1]:
            mirsetCDN2.append(edgeListCDN2[y][0])
    superMiRListCDN2.append(mirsetCDN2)
print("Second Master miRNA list complete.")

print("Creating CDN: Finding intersections of the miRNA lists")
CDN = []
for x in range(len(ceRListCDN)):
    if x == int(len(ceRListCDN)*(1/10)):
        print("10% complete...")
    if x == int(len(ceRListCDN)*(2/10)):
        print("20% complete...")
    if x == int(len(ceRListCDN)*(3/10)):
        print("30% complete...")
    if x == int(len(ceRListCDN)*(4/10)):
        print("40% complete...")
    if x == int(len(ceRListCDN)*(5/10)):
        print("50% complete...")
    if x == int(len(ceRListCDN)*(6/10)):
        print("60% complete...")
    if x == int(len(ceRListCDN)*(7/10)):
        print("70% complete...")
    if x == int(len(ceRListCDN)*(8/10)):
        print("80% complete...")
    if x == int(len(ceRListCDN)*(9/10)):
        print("90% complete...")
    list1 = superMiRListCDN1[x]
    for y in range(len(diseaseListCDN)):
        list2 = superMiRListCDN2[y]
        count = len(union(list1, list2))
        zScore = (count - meanCDN)/sigmaCDN
        if zScore > 3:
            CDN.append([ceRListCDN[x], diseaseListCDN[y], zScore])
print("CDN 100% Complete.")

# create node and edge lists for the CGN after hypergeometric reduction
print("Creating node and edge lists...")
ceRNACDNList = []
diseaseCDNList = []
for x in range(0, len(CDN)):
    if CDN[x][0] not in ceRNACDNList:
        ceRNACDNList.append(CDN[x][0])
    if CDN[x][1] not in diseaseCDNList:
        diseaseCDNList.append(CDN[x][1])
print("CDN ceRNA Nodes:", len(ceRNACDNList))
print("CDN Disease Nodes:", len(diseaseCDNList))
print("CDN Edges:", len(CDN))

# write the list as a pandas dataframe and export as a csv file
columnsCDN = ['ceRNA','Disease','Score']
dfCGN = pd.DataFrame(CDN, columns = columnsCDN)
dfCGN.to_csv('ceRNetwork_CDN.csv', index = False)

# MAKING THE MIRNA PROJECTION (CMP)

print("Creating Master miRNA Lists for Each ceRNA Node")
superMiRListCMP = []
for x in range(len(ceRNACRNList1)):
    if x == int(len(ceRNACRNList1)*0.1):
        print("10% complete...")
    if x == int(len(ceRNACRNList1)*0.2):
        print("20% complete...")
    if x == int(len(ceRNACRNList1)*0.3):
        print("30% complete...")
    if x == int(len(ceRNACRNList1)*0.4):
        print("40% complete...")
    if x == int(len(ceRNACRNList1)*0.5):
        print("50% complete...")
    if x == int(len(ceRNACRNList1)*0.6):
        print("60% complete...")
    if x == int(len(ceRNACRNList1)*0.7):
        print("70% complete...")
    if x == int(len(ceRNACRNList1)*0.8):
        print("80% complete...")
    if x == int(len(ceRNACRNList1)*0.9):
        print("90% complete...")
    tempmirlist = []
    for y in range(lenBeforeDrop1):
        if ceRNetwork1Duplicates[y] == False:
            if ceRNACRNList1[x] == ceRNetwork1['ceRNA'][y]:
                tempmirlist.append(ceRNetwork1['miRNA'][y])
    superMiRListCMP.append(tempmirlist)
print("Master miRNA List Complete.")

print("Finding Intersections of miRNA Lists for Each ceRNA-ceRNA Pair")
CMP = []
for x in range(len(superMiRListCMP)):
    list1 = superMiRListCMP[x]
    for y in range(x + 1, len(superMiRListCMP)):
        list2 = superMiRListCMP[y]
        count = len(union(list1, list2))
        if count != 0:
            CMP.append([superMiRListCMP[x], superMiRListCMP[y], count])
colsCMP = ['ceRNA1', 'ceRNA2', 'Score']
dfCMP = pd.DataFrame(CMP, columns = colsCMP)
dfCMP.to_csv('miRNA_projection.csv', index = False)

# MAKING THE GENE PROJECTION (CGP)

print("Creating Master Gene Lists for Each ceRNA Node")
superGeneListCMP = []
for item in ceRNACGNList: 
    tempgenelist = []
    for x in range(len(CGN)):
        if item == CGN['ceRNA'][x]:
            tempgenelist.append(CGN['Gene'][x])
    superGeneListCMP.append(tempgenelist)
print("Master Gene List Complete.")

print("Finding ")

# CONSTRUCTING TRIPARTITE NETWORKS: CeRNetwork (CMGN)

# create the edge lists
CMGN1 = edgeListCGN1
CMGN2 = []
for x in range(0, len(edgeListCGN2)):
    CMGN2.append([edgeListCGN2[x][1], edgeListCGN2[x][0]])
for x in range(0, len(CMGN2)):
    CMGN1.append(CMGN2[x])

columnsCMGN = ['ceRNA_mRNA', 'miRNA']
dfCMGN = pd.DataFrame(CMGN1, columns = columnsCMGN)
dfCMGN.to_csv('ceRNetwork_CMGN.csv', index = False)

# CONSTRUCTING TRIPARTITE NETWORKS: CeRDiseaseNetwork (CMDN)

# create the edge lists
CMDN1 = edgeListCDN1
CMDN2 = []
for x in range(0, len(edgeListCDN2)):
    CMDN2.append([edgeListCDN2[x][1], edgeListCDN2[x][0]])
for x in range(0, len(CMDN2)):
    CMDN1.append(CMDN2[x])

columnsCMDN = ['miRNA_Disease', 'Gene']
dfCMDN = pd.DataFrame(CMDN1, columns = columnsCMDN)
dfCMDN.to_csv('ceRNetwork_CMDN.csv', index = False)

# CONSTRUCTING TRIPARTITE NETWORKS: MGDN

# create the edge lists
MGDN1 = edgeListMDN1
MGDN2 = []
for x in range(0, len(edgeListMDN2)):
    MGDN2.append([edgeListMDN2[x][1], edgeListCDN2[x][0]])
for x in range(0, len(MGDN2)):
    MGDN1.append(MGDN2[x])

columnsMGDN = ['miRNA_Disease', 'Gene']
dfMGDN = pd.DataFrame(MGDN1, columns = columnsMGDN)
dfMGDN.to_csv('ceRNetwork_MGDN.csv', index = False)
