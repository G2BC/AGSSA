import numpy as np
import os

class UtilsTreinamento:
    def __init__(self) -> None:
        pass

    @staticmethod
    def nucleotideToCodon(base:list,nucleotides:dict,tripleNucleotides:list)->int:

        codonNumber = 0  # 0 to 63

        for p in range(3):
            if ("AGCTagct".find(tripleNucleotides[p]) != -1):
                codonNumber += base[p] * nucleotides[tripleNucleotides[p]]
            else:
                if tripleNucleotides[p] == "-":
                    codonNumber = 64  # codon com pelo menos 1
                elif tripleNucleotides[p] == "n" or tripleNucleotides[p] == "N":
                    codonNumber = 65  # codon com pelo menos 1 N
                else:
                    codonNumber = 66  # codon com caractere IUPAC
                break
        return codonNumber

    @staticmethod
    def nucleotideSequenceToCodons(base:list,nucleotides:dict,sequence:str)->list:

        tripleNucleotidesList = [sequence[i:i + 3]
                                for i in range(0, len(sequence), 3)]

        return list(map(lambda triple: UtilsTreinamento.nucleotideToCodon(base,nucleotides,triple), tripleNucleotidesList))

    @staticmethod
    def CompareResultAndAnnotation(accession, annot, clope, verbose=False):

        # REQUIREMENT: accession and annot must be compatible in size and ordering !!!

        # pegamos a lista de genotipos anotados distintos
        # <= numero de sequencias no data set de treinamento

        # definimos o nome das colunas da matriz de distribuição
        species = list(set(annot))

    #     if verbose:
    #         #print("list of genotypes in dataset:\n",species)# numeramos os genótipos segundo a ordem na lista

        # definimos o n umero de colunas
        nspecies = len(species)

        # criamos um dict com a anotacao da especie para atribuir um ordinal à cada especie anotada
        species_nmb = {}
        # atribuindo um consecutivo a cada genotipo
        ctr = 0
        for sp in species:
            species_nmb[sp] = ctr
            ctr += 1

        # calculamos o numero de transacoes de cada espécie anotada (a multiplicidade delas)
        colsum = np.zeros(nspecies, dtype=int)
        for sp in annot:
            colsum[species_nmb[sp]] += 1

    #     if verbose:
    #         #print("number of strain by genotype (colsum) \n",colsum)

        # inicializa a matriz de distribuição
        K = clope.max_cluster_number
        distr_matrix = np.zeros((K, nspecies), dtype=int)

        # constroi a matriz de distribuição
        # varre dict transaction {trans_id:cluster}
        for trans_id in range(len(clope.transaction)):  # trans_id =  clope key
            cluster = clope.transaction[trans_id]
            distr_matrix[cluster, species_nmb[annot[trans_id]]] += 1  # OK

        # attributing most probable genotype to each cluster
        annotation = {}
        for cluster in range(K):
            if np.sum(distr_matrix[cluster]) > 0:
                # anotando cada cluster com a especie mais provável
                # ORDENANDO DE MAIOR A MENOR
                sort = -np.sort(-distr_matrix[cluster])
                if sort[1] == sort[0]:  # tem mais de um máximo
                    # accession[cluster] is not good
                    annotation[cluster] = "MULTIPLE SPECIES"
                    # "UNKNOWN"  TWO OR MORE SPECIS WITH EQUAL FREQUENCY
                else:
                    annotation[cluster] = species[np.argmax(distr_matrix[cluster])]
            else:
                # CLUSTER VAZIO - RARIDADE MAS NÃO IMPOSSIVEL
                annotation[cluster] = "EMPTY CLUSTER"

        return np.array(distr_matrix), annotation

    @staticmethod
    def GetClusteringResolution(distrmat):

        # cálculo da métrica de erro relativo por cluster:
        # metric = indice_agrupamento / indice_dispersão
        # indice_agrupamento = nseqs_of_more_freq_species/nmb_of_seqs/nclust # <=1
        # indice_dispersão = Nspecies - not_present_species # >=1

        resolution = []
        n = np.sum(distrmat)  # numero total de sequencias
        dispersed_seqs = 0
        Nspecies = distrmat.shape[1]
        nclust = distrmat.shape[0]
        for cluster in range(nclust):  # varre os clusters clope
            # pega a linha da matriz de distribuição do cluster,
            row = list(distrmat[cluster])
            # que diz qtas sequencias de cada espécies estão representadas no cluster
            # conta qtas espécies não tem representante no cluster
            not_present_species = row.count(0)
            nmb_of_seqs = np.sum(row)  # number of seqs in this cluster
            # numero de seqs da espécie mais frequente no cluster
            nseqs_of_more_freq_species = np.max(row)
            # number of seqs outside the most frequented cluster
            dispersed_seqs += nmb_of_seqs - nseqs_of_more_freq_species
            if nmb_of_seqs > 0:
                indice_agrupamento = nseqs_of_more_freq_species/nmb_of_seqs/nclust**0.5  # <=1
                indice_dispersão = Nspecies - not_present_species  # >=1
                resolution.append(indice_agrupamento / indice_dispersão)
            else:
                resolution.append(0)  # EMPTY CLUSTER CASE ------------->

        return np.array(resolution), 1 - dispersed_seqs/n

    @staticmethod
    def get_accession(new_accession, transactions):
        for j in range(len(transactions)):
            if new_accession in transactions[j][0].replace("/", ""):
                return transactions[j][0]
        return "unknown"
    
    @staticmethod
    def MeanResolution(x, w):
        WM = []
        n = len(x)
        if len(w) != n:
            # print("WRONG WEIGHTS ----------------------------------------------")
            return -1
        m = np.sum(w)
        for i in range(n):
            WM.append(x[i] * w[i]/m)
        return WM


    # calcula a distancia entre 2 vetores categóricos:
    @staticmethod
    def categoricaldist(u, v, var, mode=0):
    # numero de entradas distintas / total de entradas
        dist = 0
        for i in range(len(u)):
            if u[i] != v[i]:
                if mode == 0:
                    # weighted -> INVERSAMENTE PROPORCIONAL À VARIABILIDADE DO SITE
                    dist += 1/var[i]
                else:
                    dist += 1  # SEM CONSIDERAR A VARIABILIDADE DO SITE

        return dist/len(u)
    


    # Função para converter um ClusterNode em Newick
    @staticmethod
    def get_newick(node, labels):
        if node.is_leaf():
            return labels[node.id]
        left = UtilsTreinamento.get_newick(node.left, labels)
        right = UtilsTreinamento.get_newick(node.right, labels)
        return f"({left},{right})"

    # Função personalizada para obter a representação Newick com comprimento dos ramos
    @staticmethod
    def get_newick_with_branch_length(node, lista_de_nomes):
        if node.is_leaf():
            return f"{lista_de_nomes[node.id]}:{node.dist}"
        else:
            esquerda = UtilsTreinamento.get_newick_with_branch_length(node.left, lista_de_nomes)
            direita = UtilsTreinamento.get_newick_with_branch_length(node.right, lista_de_nomes)
            return f"({esquerda},{direita}):{node.dist}"
        
    @staticmethod
    def zip_directory(directory, zip_filepath):
        try:
            base_dir = os.path.basename(directory)
            parent_dir = os.path.dirname(directory)
            os.system(
                f"cd '{parent_dir}' && zip -r '{zip_filepath}.zip' '{base_dir}'")
        except Exception as e:
            print(f'Erro ao compactar diretório: {e}')