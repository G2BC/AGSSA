# import our three important packages
#!/usr/bin/env python
# -*- coding: utf-8 -*-
from collections import Counter
import os
import sys
import shutil
import Bio as Bio
from Bio import SeqIO
from Bio import AlignIO

import time as time
import numpy as np

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform

import pandas as pd
import random
from io import StringIO

import dill  # for model object handling
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio.Phylo.BaseTree import Clade as Clade
from Bio.Phylo.BaseTree import Tree as BioTree
from collections import defaultdict
from .Untils_treinamento import UtilsTreinamento 
from .Model import Model 
from .clusterizer import CLOPE, Cluster



# from ete3 import Tree, TreeStyle

start_time_main = time.time()

# Caminho para as sequências
fname = sys.argv[1]
# Caminho do arquivo de anotações
arquivo_anotacoes = sys.argv[2]
# Caminho para os arquivos de saída
output_path = sys.argv[3]
# Nome da espécie
especie = sys.argv[4]


# Nome da espécie
os.makedirs(output_path,exist_ok=True)

especie ='spike'

# FASE 1 - PREPROCESSAMENTO
# # LEITURA DO DATASET DE TREINAMENTO & IDENTIFICAÇÃO/CODIFICAÇÃO DE CLASSES PRIMÁRIAS

base = [16, 4, 1]

nucleotides = {
    "A": 0,
    "G": 1,
    "C": 2,
    "T": 3,
    "a": 0,
    "g": 1,
    "c": 2,
    "t": 3,
}

stops = ["TAA", "TAG", "TGA"]

stopnmb = []

for triplet in stops:
    codonNumber = 0  # 0 to 63
    for p in range(3):
        codonNumber += base[p] * nucleotides[triplet[p]]
    stopnmb.append(codonNumber)

# print("stop codon ordinals ", stopnmb)

# import the sequences we will use. These are 16s sequences from GenBank
# example: https://www.ncbi.nlm.nih.gov/nuccore/FJ039971.1?report=genbank

t = []
with open(fname) as handle:
    for values in SeqIO.FastaIO.SimpleFastaParser(handle):
        t.append(values)
        # print(t[-1][0])

Nseq = len(t)

# print("number of sequences in file ", fname, " = ", Nseq)

# Open the alignment file AGAIN BUT NOW as a MultipleSeqAlignment object

with open(fname, "r") as aln:
    alignment = AlignIO.read(aln, "fasta")

idlst = []
sequences = []
labels = []
ctr = 0  # conta não redundantes
ctrdup = 0  # conta sequencias duplicadas
for strain in alignment:
    #     #print(strain.description)
    if 'accn|' in strain.description:
        accn = strain.description.split("accn|")[1]
    else:
        accn = strain.description

    if accn not in idlst:
        idlst.append(accn)
        # strain.description # o id só saia cortado e acahava depois que itnha sequncias repetidas
        labels.append("strain"+str(ctr))
        sequences.append(alignment[ctr].seq)
        ctr += 1
    else:
        # print("WARNING: duplicated description ", accn)
        ctrdup += 1

Nseqs = len(sequences)
# print("number of non-redundant accessions in the alignment", Nseqs, " (# replicated: ", ctrdup, ")")

# transformando sequencias de nucleotideos a codons


sequenceLength = len(sequences[0])

if (sequenceLength % 3 != 0):
    # print(sequenceLength)
    # print("ERRO: tamanho não Múltiplo de 3")
    exit(1)
elif (next((True for seq in sequences if len(seq) > sequenceLength), False)):
    pass
    # print("ERRO: as sequências não tem o mesmo tamanho")
else:
    nprm = int(sequenceLength / 3)  # comprimento das seqs em codons

# print(f"O dataset contem {Nseqs} sequências com {nprm} resíduos.")
# print("Transformando as sequências de nucleotídeos em sequências númericas de códons... ")

start = time.time()




codonSequences = [UtilsTreinamento.nucleotideSequenceToCodons(base = base, nucleotides = nucleotides, sequence = seq) for seq in sequences]



# print(f"Transformação concluída em {(time.time() - start) / 60} minutos")

# check codon numbering
# print("minimum codon number: ", np.min(codonSequences), "maximum codon number: ", np.max(codonSequences))

# filtrando sequências com stop codons TGA,TAG,TAA e tirando do final
filtered_codonSeqs = []
filtered_Seqs = []
filtered_SeqID = []
filtered_Label = []
ctr = 0
laststop = 0

#lucas alterando
#Guarda id das sequencias excluidas para atualizar anotações depois
seq_excluded = []


for seq in codonSequences:  # olho: se assume que todas as sequencias tem ou não tem stopcodon no final
    # se não for assim, o corte não será bem feito e teremos sequencias de comprimento diferente
    # podendo dar erro mais na frente
    stopfound = False
    for stp in stopnmb:
        stoppos = np.where(np.array(seq) == stp)[0]
        if len(stoppos) > 0:
            if len(stoppos) == 1 and (stoppos[0] == nprm-1):  # last stop found
                laststop = 1
            else:
                print(" stop codon ", stp, " found in frame in sequence ", alignment[ctr].id, " at positions ", stoppos[0])
                seq_excluded.append(alignment[ctr].id.strip())
                stopfound = True
                break

    if not stopfound:  # somente append as que nao tem stop codons in frame
        filtered_SeqID.append(idlst[ctr])  # new 05.10.2023
        filtered_Label.append(labels[ctr])
        filtered_Seqs.append(sequences[ctr])
        if laststop > 0:
            # excluindo o ultimo stop se existir
            filtered_codonSeqs.append(seq[:-laststop])
        else:
            filtered_codonSeqs.append(seq)
    ctr += 1




# check stop codon filtering

Nseq = len(filtered_codonSeqs)
# print("number of sequences before filtering: ", len(codonSequences), " after: ", Nseq)
# print("length in codons of sequences before filtering: ", len(codonSequences[0]), " after: ", len(filtered_codonSeqs[0]))

# transformando o dataset em matriz de codons

cdnseq = np.array(filtered_codonSeqs)
# print(cdnseq.shape)

# contando número de codons distintos por coluna (sítios)
cdns = []
V = []

for pos in range(cdnseq.shape[1]):
    cdnlst = list(set(list(cdnseq[:, pos])))
    V.append(len(cdnlst))
    cdns.append(cdnlst)

V = np.array(V)

# ordeno as colunas de maior a menor variabilidade
# Idx=np.argsort(-V)
# #print("Idx",Idx)
# #print("V[Idx]",V[Idx])

# seleciono as colunas com variabilidade > 1
var_sites = np.argwhere(V > 1)
Var_Sites = []
for site in var_sites:
    Var_Sites.append(int(site))

var = V[var_sites]

# #print("var",var)

IDsz = len(var_sites)
# print("ID size ", IDsz, f" % of informative codons: {100*IDsz/nprm:.2f}")


# gero os IDs com os codons de cada seq (lista de listas)

ID = cdnseq[:, var_sites].tolist()  # vetor coluna de IDs como vetores linha

strID = []  # lista de IDs como strings
for lst in ID:
    strID.append(str(lst).replace('[', '').replace(']', ''))

# CBUC
# gero a lista de string IDs das classes primárias - unique

# gero a lista de string IDs das classes primárias - unique
start = time.time()
p = []  # vetor com string das IDs das classes primárias
pvec = []
pClass = []
pAcc = []
pSeq = []
P = 0
ctr = 0
for Id in strID:
    if Id not in p:  # filtro sequencias iguais
        p.append(Id)
        pvec.append(np.array(Id.split(", ")))
        pClass.append(P)  # numero a classe primaria
        pAcc.append(filtered_SeqID[ctr])  # capturo o accession
        pSeq.append(filtered_Seqs[ctr])  # capturo o accession
        P += 1
    ctr += 1

# print("check primary class 0 and its accession, string ID and char vector ID and FASTA")
# print("NUMBER OF THE PRIMARY CLASS: ", pClass[0], " ACCESSION: ", pAcc[0])
# print("STRING ID: \n", p[0])
# print("VECTOR ID: \n", pvec[0])
# print("NUCLEOTIDE SEQUENCE: \n", pSeq[0])

# numero de classes primarias
# print("there were found ", P, " primary classes in the dataset with ", Nseqs, " sequences in ", time.time()-start, " segundos")
# print(f"different sequences: {P/Nseq*100:.2f} %")

# check
# print(len(pvec), "=", P)

# criando dataframe com o Gene Variant File

UniqueSeqs = pSeq
UniqueSeqIDs = pAcc
UniqueClass = pClass

gvf = pd.DataFrame({'class': pClass, 'accession': pAcc,
                   'fasta_aln': pd.Series(UniqueSeqs)})  # ,index=0)

gvf

# controle de qualidade do dataset de treinamento
# quantas sequencias tem de cada classe primária?
# >>>>> deve ser 1 só

# histograma de frequencia das classes primárias - depende do ordenamento

counts = Counter(pClass)
# #print(counts)
numero_mais_comum, frequencia_mais_comum = counts.most_common(1)[0]

if frequencia_mais_comum > 1:
    # print("WARNING: THERE ARE REDUNDANT SEQUENCES IN THE TRAINING DATASET ")
    # Exiba o número mais comum e sua frequência
    # print(f"The most common class is {numero_mais_comum} with {frequencia_mais_comum} copies")

    plt.figure(figsize=(14, 2))
    plt.hist(pClass, P)
    plt.grid()
    plt.xlabel("primary class")
    plt.ylabel("# of seqs by primary class")
    # plt.show()

else:
    # print("NO REDUNDANT SEQUENCES HAVE BEEN FOUND IN THE TRAINING DATASET - PASSED ...")
    pass



model = Model(especie, nmbofclasses=P)
model.ListOfVarSites = Var_Sites  # CBUC

# FASE 2 - ANÁLISE QUANTITATIVA

# 1. AGRUPAMENTO
# 2. ANOTAÇÃO DAS CLASSES PRIMÁRIAS
# 3. CÁLCULO DA RESOLUÇÃO DO CLASSIFICADOR


# acessando à anotação de cada classe primária


# Dicionário para armazenar as anotações
annot = {}

# Ler o arquivo de anotações
with open(arquivo_anotacoes, 'r') as anotacoes_file:
    for linha in anotacoes_file:
        
        id_seq_match = False
        for seq in seq_excluded:
            if seq in linha:
                id_seq_match = True
                break
        if id_seq_match:
            continue
        # Dividir a linha usando a vírgula como separador
        partes = linha.strip().split(',')  # primeira parte ID longo terminando en ACCN,
        # segunda parte o tipo atribuido pelo anotador externo
        # Verificar se a linha tem o formato correto
        if len(partes) == 2:
            cabecalho, anotacao = partes
            annot['_'.join(cabecalho.split('_')[1:])
                  ] = anotacao  # new 20.11.2023
            # annot[cabecalho.split("_")[:1]] = anotacao # new 05.10.2023 / 20/10/2023 Alterado o formato do arquivo de anotações
           
        else:
            print(f"A linha '{linha}' não possui o formato esperado.")
            

# Exibir o dicionário de anotações
nAnnotSeq = len(annot)
if nAnnotSeq != Nseq:
    print(
        f"PROBLEM: THE NUMBER OF FILTERED SEQUENCES({nAnnotSeq}) AND ANNOTATED SEQUENCES({Nseq}) DOES NOT MATCH ...CALL DIEGO!")
else:
    # print("we read ", nAnnotSeq, " annotated sequence accessions\n", annot)  # {accesion: Pangolin genotype}
    pass

# check - all input sequence must be annotated

# IDENTIFYING DIFFERENT ANNOTATED SPECIES (SPECIES)

# pegando as espécies virais anotadas
species = list(set(annot.values()))
Nspecies = len(species)
# print("we found ", Nspecies, " annotated genotypes\n", species)

# funções que constroem a matriz de distribuição e calculam a resolução do agrupamento usando o ground truth (supervisão para ajuste do hyper-parâmetnro "repulsion")

# funções que constroem a matriz de distribuição e calculam o erro usando o ground truth

# CompareResultAndAnnotation retona a "distr_matrix" e "annotation"
# OS CLUSTERS SÃO ANOTADOS COM 3 CLASSES:
# "EMPTY CLUSTER", quando está vazio,
# "MULTIPLE SPECIES" quando 2 ou mais espécies tem o mesmo numero de sequencias no cluster e não
# é possível definir a espécie (pre)dominante no cluster, e
# "NOME DA ESPECIE DOMINANTE", quando uma única espécie é a mais frequente no cluster
# (o que não significa que seja a única...)



# OPTIMIZATION PROCESS OF THE UNSUPERVISED
# CLUSTERING ALGORITHM HYPER-PARAMETER
# it is done by:
# 1. varying repulsion in a given search interval
# 2. doing several random iterations of CLOPE clustering for each repulsion sorting input data
# 3. recording the clustering that maximizes our quality criteria which is based on less
# dispersion of annotated especies among the found clusters
# >>>> A quality threshold is defined as an early stopping trigger

# MAIN -----------------------------------------------------


# pegando a lista de accession
Acc = gvf['accession']
# tirando barras para match com annot
Accession0 = []
for acc in Acc:
    Accession0.append(acc.replace("/", ""))  # um dataset de Inês tinha "/"
# #print(Accession0[:3])

# pegando o groundthruth da anotação
GroundTruth0 = list(annot.values())
# GroundTruth0[:3]

# pegando os vetores de chars dos códigos das classes primárias como transações
InputData0 = pvec
nmb_trans = len(pvec)

nmb_annot = len(GroundTruth0)
# print('check: number of transactions = annotations: ', nmb_trans, '=', nmb_annot)

if nmb_trans != nmb_annot:
    # print("DATASET NOT FULLY OR WRONGLY ANNOTATED .....\nRESULTS MAY BE CONFUSING....")
    pass

noiseLimit = 0.0  # dummy

nRandIterations = 30  # 50 equivalent to bootstrap
ClusterResolutionThres = 0.5  # 1.00

# print("Looking for optimal repulsion with ", nRandIterations, "random iterations and threshold of minimum cluster resolution of ", ClusterResolutionThres)

optimum = False
minRep = 1
maxRep = 16
stpRep = 2

trans_len = len(InputData0[0])

best_repulsion = -1
best_resolution = 0
K = nmb_trans

# BUSCANDO PELA MELHOR REPULSION ACOMPANHANDO NOSSA METRICA DE AGRUPAMENTO

tic = time.time()

# >= 11 da 1005 acc de clusters com 200 iteracoes randomicas
for repulsion in range(minRep, maxRep, stpRep):

    # print("\n>>>> repulsion", repulsion)

    cluster_resolution = []
    Resolution = []

    best_K = nmb_trans

    codes = InputData0
    GT = GroundTruth0

    # EXECUTANDO nRandIterations AGRUPAMENTOS COM CLOPE (NÃO SUPERVIDIONADO)
    # COM O VALOR DE REPULSION DADO, OLHANDO SEMPRE PARA O MELHOR AGRUPAMENTO

    for i in range(nRandIterations):

        seed = None  # 150 <<<<<<<<<<<<<<<<<<<<<<<< ORDENAÇÂO ALEATORIA
        np.random.seed(seed)

    # not needed because clope shuffle input data every it is called << NEW 07.10.2023
    #     InputData=InputData0.copy()
#         list_zipped = list(zip(Accession0,InputData0,GroundTruth0))
#         np.random.shuffle(list_zipped)
#         Accession, InputData, GroundTruth = zip(*list_zipped) #unzipping

        # transferimos os dados originais para as vari[aveis de trabalho

        Accession = Accession0.copy()
        InputData = InputData0.copy()
        GroundTruth = GroundTruth0.copy()

        # #print('first transactions:\n',InputData[:3])
        # #print('first annotations:\n',GroundTruth[:3])

        # criamos o dataset formato CLOPE ordeando de forma aleatoria os dados de entrada

        Data = {}  # usamos um dicionario para os dados com as features indexadas por posição

        for trans_nmb in range(0, len(InputData)):
            # inicializa a transação efetiva
            Data[trans_nmb] = [''] * trans_len
            fst = 0  # for SARS2
            missing = []
            for index in range(fst, len(InputData[trans_nmb])):
                # adicionando a posição (coluna) da feature se não for missing
                if InputData[trans_nmb][index] != '?':
                    Data[trans_nmb][index] = str(
                        index) + ":"+InputData[trans_nmb][index].replace(" ", "")  # add feature position
                else:  # counting missings
                    missing.append(index)
                if len(missing) > 0:
                    # print("WARNING: there were ", len(missing), " missing features in transaction ", trans_nmb, " at positions:\n", missing)
                    pass

        # criamos um novo clusterizador
        clope = CLOPE(print_step=1000, random_seed=seed)

        # Dados iniciais # Inicializamos o algoritmo
        clope.init_clusters(Data, repulsion, noiseLimit)

        # Iteramos até o metodo de clusterização não supervisionado convergir
        ctr = 0
        while clope.next_step(Data, repulsion, noiseLimit) > 0:
            ctr += 1
            clope.print_history_count(repulsion, seed)

        print("RANDOM ITERATION ", i+1, " -> clope iterations: ",
              ctr, " number of clusters: ", clope.max_cluster_number)

#         #print("CLOPE number of clusters ",clope.max_cluster_number)
        # #print(clope.transaction.values())
        # #print("annot",list(GroundTruth))

        # PARTE SUPERVISIONADA ================================================================

        # construção da matriz de distribuição: numero de cada genotipo (colunas) em cada cluster clope
        # (linhas)

        distrmat, map_species2cluster = UtilsTreinamento.CompareResultAndAnnotation(
            Accession, GroundTruth, clope, True)

        # cálculo da "resolução da classificação" ( novo conceito nosso: a resolução aumenta
        # na medida que os clusters são mais homogêneos no sentido de conter um mesmo genótipo
        # preferencialmente )
        # A maior homogeneidade é qdo todas as sequencias de um cluster clope pertencem ao mesmo
        # genotipo. não importa que vários clusters clope contenha esse mesmo genótipo. o que conta
        # é a homegeneidade media * homgeneidade minima dos clusters

        resolution_lst, tot_resolution = UtilsTreinamento.GetClusteringResolution(distrmat)

        cluster_resolution.append(resolution_lst)
        Resolution.append(tot_resolution)

        # condition to choose the best clustering during iterations USING MIN*MEAN Resolution
        resol = np.sqrt(
            np.min(cluster_resolution[-1])*np.mean(cluster_resolution[-1]))

        # GUARDA AS INFORMAÇÕES DA MELHOR CLUSTERIZAÇÃO ATÉ O MOMENTO
        if (resol > best_resolution):
            # and (len(distrmat) < best_K): CRITERIO QUE USEI SEM MUITO CRITÉRIO , KKKKKKKK
            # A CONDIÇÃO ACIMA BUSCA A COMBINAÇÃO COM MENOR DISTRIBUIÇÃO (MAIOR ACURACIA)
            # E MENOR NUMERO DE CLUSTERS !!!!!! -->>>>>> tirei essa condição 18.10.2023
            best_clusters_resolution = cluster_resolution[-1]
            best_resolution = resol  # or mean USING MIN Resolution
            best_K = len(distrmat)
            K = best_K
            best_distrmat = distrmat
            best_annotation = map_species2cluster
            map_class2cluster = clope.transaction
            codes = InputData
            GT = GroundTruth
            best_accession = Accession
            best_repulsion = repulsion
            # print(">>>>>>>>>>>>>>>> best found with ", best_K, " clusters. resolution = ", best_resolution)

        del clope  # deleta o classificador - economia de memoria

    # controle de optimalidade
    # AQUI SE DEFINE QUÃO EXIGENTE É O ESTUDO, SE PERMITE ALGUM NIVEL DE DISTRIBUIÇÃO OU NÃO
    # SE  ClusterResolutionThres = 1 => É PORQUE NÃO HOUVE NENHUMA MSITURA DE ESPECIES NOS CLUSTERS, OU SEJA
    # CADA CLUSTER CONTINHA UM UNICO TIPO DE ESPECIE ANOTADA
    if best_resolution >= ClusterResolutionThres:  # NEW 08.10.2023
        # print("satisfied optimality criteria")
        optimum = True
        break  # saimos do loop aumentando repulsion

# print("OPTIMIZATION ENDED\nNumber of repulsions tested: ", int((maxRep-minRep)/stpRep), "\nNumber of random iterations: ", nRandIterations, "\nElapsed time: ", (time.time()-tic)/60, " min")

# PREENCHENDO O MODELO COM A CLUSTERIZAÇÃO ÓTIMA (MELHOR)


# -----------------------------------------------------------------------------------------------------


if optimum:
    # print("optimum criterium found for repulsion = ", best_repulsion, " (min*mean)^(1/2) resolution = ", best_resolution)
    pass
else:
    # print("optimum not found for repulsion in the interval [", minRep, ",", maxRep, "] best (min*mean)^(1/2) resolution = ", best_resolution)
    pass

# evaluate clustering resolution
x = best_clusters_resolution
w = np.sum(best_distrmat, axis=1)

# print("Check", len(x), len(w))

wr = UtilsTreinamento.MeanResolution(x, w)

model.NmbOfClusters = K
model.Resolution = wr  # mean
model.MinClusterResolution = np.min(wr)  # MODEL
model.DistributionMatrix = best_distrmat  # MODEL
model.CodeOfClass = codes  # origin CBUC but randomly shuffled by AGUA
model.GroundTruth = GT  # ANNOTATION
model.TheClusterOfClass = map_class2cluster  # CLOPE
model.TheSpeciesOfCluster = list(
    best_annotation.values())  # CLOPE + ANNOTATION
model.CLOPE_repulsion = best_repulsion  # MODEL
baccession = []
for acc in best_accession:
    baccession.append(UtilsTreinamento.get_accession(acc, t))
model.Accession = baccession  # CBUC but shuffled by AGUA

# print("Mean resolution of the model ", np.sum(model.Resolution), " standard deviation of resolution ", np.std(model.Resolution))
# print("Number of Classes :", model.NmbOfClasses)
#     #print("CodesOf Classes :\n",model.CodeOfClass)
#     #print("Annotation of Classes :\n",model.GroundTruth)
# print("Number of clusters :", model.NmbOfClusters)
df_distrmat_SpeciesByClusters = pd.DataFrame(
    zip(model.DistributionMatrix, model.TheSpeciesOfCluster))
# print("Confusion matrix\n", df_distrmat_SpeciesByClusters)

# visualizing clusters with multiple species

# visualizing clusters with multiple species

# ESSES CLUTERS DEVEM CONTER SEQUENCIAS QUE NÃO FORAM BEM CLASSIFICADAS ...

mixed_clusters_nmb = 0
for i in range(len(df_distrmat_SpeciesByClusters[1])):
    if df_distrmat_SpeciesByClusters[1][i] == "MULTIPLE SPECIES":
        mixed_clusters_nmb += 1
        # print(i, df_distrmat_SpeciesByClusters[0][i], " nseqs (", np.sum(df_distrmat_SpeciesByClusters[0][i]), ")")

# FASE2 - CÁLCULO DA MATRIZ DE DISTÂNCIAS ENTRE PARES DE CLASSES PRIMÁRIAS

# calculating distances between pairs of IDs -> put in categorical distance matrix





# -------------------------------------------------------------------


# print("FASE2 - CÁLCULO DA MATRIZ DE DISTÂNCIAS ENTRE PARES DE CLASSES PRIMÁRIAS ")

fvec = model.CodeOfClass

start = time.time()

cat_distmatrix = []
for i in range(P-1):
    cat_distmatrix.append([])
    for j in range(i+1, P):
        #         dotdistmatrix[-1].append(np.dot(pvec[i],pvec[j]))
        # calling weighted
        cat_distmatrix[-1].append(UtilsTreinamento.categoricaldist(fvec[i], fvec[j], var, 0))
        # lçast prm above: calling in mode=0 sem pesos, calling in mode=1 com pesos
#         #print(i,j,np.dot(pvec[i],pvec[j]))

# check
# for distvec in cat_distmatrix:
#     #print(distvec)

# gerando matriz de distancias simetrica

# nota: os vetores DistVec , fst e scd estão alinhados, ou seja DistVec[k] é a distância entre a classe fst[k] e a scd[k] !!!

DistMat = np.zeros((P, P), dtype=float)
labels = []
# fst & scd são o par de indexadores de cada distância (as classes primárias que são comparadas)
# atualizando labels com a especie predita
fst = []
scd = []
for i in range(P-1):
    labels.append("pclass:"+str(i) +
                  "|pacc:"+pAcc[i] +
                  "|cluster:"+str(model.TheClusterOfClass[i]) +
                  "|annot:"+model.GroundTruth[i] +
                  "|pred:"+model.TheSpeciesOfCluster[model.TheClusterOfClass[i]])
    for j in range(i+1, P):
        fst.append(i)
        scd.append(j)
        DistMat[i][j] = cat_distmatrix[i][j-i-1]
        DistMat[j][i] = DistMat[i][j]  # building a symmetric matrix

# labeling last row
labels.append("pclass:"+str(P-1) +
              "|pacc:"+pAcc[P-1] +
              "|cluster:"+str(model.TheClusterOfClass[P-1]) +
              "|annot:"+model.GroundTruth[P-1] +
              "|pred:"+model.TheSpeciesOfCluster[model.TheClusterOfClass[P-1]])

# print("matriz de distâncias\n", DistMat)

# #print(DistMat.shape)

# CONSTRUINDO MATRIZ DE DISTANCIAS LOW TRIANGULAR E VETOR DE LABELS PARA GERAR ARVORE NEWICK
LowTriang = []
LowLabels = []
matsz = len(DistMat)
for i in range(matsz-1, -2, -1):
    row = []
    for j in range(i+1, matsz):
        #         #print(i,j)
        row.append(DistMat[i, j])
    if len(row) > 0:
        LowTriang.append(row)
        LowLabels.append(pAcc[i])

# print("number of samples in the distance matrix: ", len(LowLabels), "\nLabels: ", LowLabels)

# END LOW TRIANGULAR WORK

# SALVA A MATRIZ DE DISTANCIAS
np.savetxt(f'{output_path}/{model.name}_dist_mat.txt', DistMat)

# saida para a construcao de arvores com dendogram
DistVec = squareform(DistMat)
# print("distance vector:\n", DistVec)

# print(f"tempo:{time.time()-start} segundos")

# FASE 3 - ANÁLISE QUALITATIVA+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Povoando a variável contador com os dados necessários
# Cluster, Tipos, Quantidade

contador = [{} for _ in range(model.NmbOfClusters)]

# Iterar sobre cada cluster
for target in range(model.NmbOfClusters):
    # Iterar sobre cada pclass no cluster atual
    for pclass in model.TheClusterOfClass.keys():
        # Verificar se a pclass pertence ao cluster atual
        if model.TheClusterOfClass[pclass] == target:
            # Obter o valor correspondente de model.GroundTruth[pclass]
            valor = model.GroundTruth[pclass]
            # Atualizar a contagem para o valor no dicionário do cluster atual
            if valor in contador[target]:
                contador[target][valor] += 1
            else:
                contador[target][valor] = 1


# vendo a composição detalhada de cada cluster 'clope'
# #print(model.TheClusterOfClass)+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Label4Tree = {}
# for target in range(model.NmbOfClusters):
    # print(">>>>>>>> cluster ", target, " contains sequences of the following ",
  #       np.sum(model.DistributionMatrix[target]), " genotypes:")
   # for pclass in model.TheClusterOfClass.keys():
    #   if model.TheClusterOfClass[pclass] == target:
    #       #             indice = pClass.index(pclass)
            #             accession_correspondente = pAcc[indice]
       #     #print("primary class ", pclass, " with accession ", pAcc[pClass.index(pclass)],
        #          " annotated as ", model.GroundTruth[pclass])
       #     Label4Tree[pAcc[pClass.index(pclass)]] = pAcc[pClass.index(
       #         pclass)]+":"+str(target)+":"+model.GroundTruth[pclass]

# Salva em um arquivo csv os dados montados na célula anterior

# Cria cabeçalho do arquivo csv
mensagem = "cluster;qtd_seqs;tipos\n"

# Percorre o contador com a quantidade de sequências por tipo dentro de cada cluster
for i, contagem in enumerate(contador):
    mensagem += f"{i};{np.sum(model.DistributionMatrix[i])};"
    tipos = ''
    for chave, valor in contagem.items():
        tipos += f"{chave},"
    mensagem += tipos[:-1]
    mensagem += "\n"

# Gera um DataFrame com a mensagem montada
df = pd.read_csv(StringIO(mensagem), sep=';')

# Salva o arquivo csv com as informações
df.to_csv(f'{output_path}/tipos_por_clursters.csv', index=False, sep=';')

# for i, contagem in enumerate(contador):
#     mensagem = f"Cluster {i}: Total Seqs: {np.sum(model.DistributionMatrix[i])}"
#     for chave, valor in contagem.items():
#         mensagem += f" Tipo: {chave} - Qtd Sequências: {valor},"
#     #print(mensagem[:-1])


# vendo a composição detalhada de cada cluster 'clope'
# #print(model.TheClusterOfClass)
Label4Tree = {}
for target in range(model.NmbOfClusters):
    # print(">>>>>>>> cluster ",target," contains sequences of the following ",\np.sum(model.DistributionMatrix[target]), " genotypes:")
    for pclass in model.TheClusterOfClass.keys():
        if model.TheClusterOfClass[pclass] == target:
            #             indice = pClass.index(pclass)
            #             accession_correspondente = pAcc[indice]
            # print("primary class ",pclass," with accession ",pAcc[pClass.index(pclass)],\" annotated as ", model.GroundTruth[pclass])
            Label4Tree[pAcc[pClass.index(pclass)]] = pAcc[pClass.index(
                pclass)]+":"+str(target)+":"+model.GroundTruth[pclass]

# TEMPO DE PROCESSAMENTO

end_time_main = time.time()
# print(f"Treinamento Realizado em {end_time_main - start_time_main} segundos")
print(
    f"Treinamento Realizado em {(end_time_main - start_time_main) / 60} minutos")
# print(f"Foram utilizadas {Nseqs} sequências")

# Label4Tree

# visualizando árvore filogenética

# linkage methods
# The following are methods for calculating the distance between the newly formed cluster  and each .

# 1. method=’single’  also known as the Nearest Point Algorithm.

# 2. method=’complete’  also known by the Farthest Point Algorithm or Voor Hees Algorithm.

# 3. method=’average’ also called the UPGMA algorithm.

# 4. method=’weighted’ also called WPGMA.

# 5. method=’centroid’  also known as the UPGMC algorithm.

# 6. method=’median’  also known as the WPGMC algorithm.

# 7. method=’ward’ also known as the incremental algorithm.

# methods=['single', 'complete','average','weighted','centroid','median','ward']
# methodknownas=['NearestPoint', 'FarthestPoint','UPGMA','WPGMA','UPGMC','WPGMC','Incremental']
methods = ['average']
methodknownas = ['UPGMA']
# methods=['ward','average','weighted']
# methodknownas=['Incremental','UPGMA','WPGMA']
nmeth = len(methods)


# construindo arvores com distintos métodos de linkagem da scipy.cluster.hierarch.linkage
# mostrando onde os accessions foram parar sem analisar nem comparar
# para usar com tree_analyzer !!!!
ctr = 0

working_labels = labels.copy()
short_labels = [label.split('|pacc:')[1] for label in working_labels]
short_labels = [label.split('|')[0] for label in short_labels]

# #print("old: \n",short_labels)

nshort_labels = []
for i in range(len(short_labels)):
    nshort_labels.append(Label4Tree[short_labels[i]])

for meth in methods:
    linkage_matrix = linkage(DistVec, meth)

    fig, ax = plt.subplots(1, 1, figsize=(10, 18))
    fig.set_size_inches(10, 18, forward=True)

    dn = dendrogram(linkage_matrix, ax=ax, labels=nshort_labels,
                    show_leaf_counts=True, show_contracted=True,
                    orientation='left', leaf_font_size=5)
    plt.title(f"{especie} - Host independent study with linkage method '" +
              meth+"' ("+methodknownas[ctr]+")")

    plt.tight_layout(pad=4)
    # plt.show()
    tree_name = 'tree_'+methodknownas[ctr]
    fig.savefig(f'{output_path}/{tree_name}.png', dpi=fig.dpi)

    ctr += 1


# MAURICIO - AQUI TEM QUE DAR UM JEITO DE SALVAR O DENDOGRAMA OU A LINKAGE MATRIX EM FORMATO NEWICK

    # E DEPOIS LER O ARQUIVO E PLOTAR PAR CONFERIR QUE FUNCIONA

#     # Converter a matriz de ligação em uma árvore hierárquica
#     tree = to_tree(linkage_matrix, rd=False)

#     # Salvar a árvore no formato Newick
#     with open(tree_name + '.newick', 'w') as f:
#         newick = UtilsTreinamento.get_newick(tree, nshort_labels)
#         f.write(newick + ";")

#     # Carregar a árvore a partir do arquivo Newick
#     tree = Tree(tree_name + '.newick', quoted_node_names=True)

#     # Criar um estilo para a plotagem da árvore (opcional)
#     ts = TreeStyle()
#     ts.show_leaf_name = True
#     ts.show_branch_length = True

#     # Plotar a árvore
#     tree.render(tree_name + ".png", tree_style=ts)

# print("END")

# CRIANDO E EXPORTANDO ARVORE NEWICK USANDO UPGMA

# MAURICIO - AQUI TENTEI DESSA FORMA MAS VEJA QUE A ARVORE ESTA CORROMPIDA, NAO SERVE

# CRIANDO E EXPORTANDO ARVORE NEWICK USANDO UPGMA

plt.rc('font', size=10)  # Substitua 12 pelo tamanho desejado

nLowLabels = []
for i in range(len(LowLabels)):
    nLowLabels.append(Label4Tree[LowLabels[i]])

# Crie uma matriz de distâncias Biopython
dm = DistanceMatrix(nLowLabels, LowTriang)

# Use o construtor de árvore de distância UPGMA
constructor = DistanceTreeConstructor()
# tree_NJ = constructor.nj(dm)  # Substitua por 'nj' para usar o método Neighbor Joining
# Substitua por 'nj' para usar o método Neighbor Joining
tree_UPGMA = constructor.upgma(dm)

# Salve a árvore em formato Newick
# with open('arvore_NJ.nwk', 'w', encoding='utf-8') as f:
#     Phylo.write(tree_NJ, f, 'newick')

with open(f'{output_path}/arvore_UPGMA.nwk', 'w', encoding='utf-8') as f:
    Phylo.write(tree_UPGMA, f, 'newick')
# Plote a árvore (opcional)
# Crie uma figura com tamanho personalizado

# Crie uma figura com tamanho personalizado
fig, ax = plt.subplots(figsize=(18, 55))  # Ajuste o tamanho conforme desejado

# Plote a árvore
# Phylo.draw(tree_NJ, axes=ax, label_func=lambda x: '' if not x.is_terminal() else x.name)
# Phylo.draw(f'{output_path}/{tree_UPGMA}', axes=ax,
#            label_func=lambda x: '' if not x.is_terminal() else x.name)
# Exiba o gráfico
# plt.show()


# UtilsTreinamento.print_model(model, 2)

# CREATING THE GENE VARIANT FILE
# TO DO

# saving and loading the  model

# Salva o arquivo
model_filename = f"model.{model.name}.obj"
model_filepath = os.path.join(output_path, model_filename)

with open(model_filepath, "wb") as file:
    dill.dump(model, file)

# Recarregar o arquivo
with open(model_filepath, "rb") as file:
    loaded_model = dill.load(file)



# Zipando os arquivos gerados em output_path


zip_filename = 'resultados'
zip_filepath = os.path.join(output_path, zip_filename)

UtilsTreinamento.zip_directory(output_path, zip_filepath)
