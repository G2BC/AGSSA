# Sequência de passos para analise com o AGUA:
# 1. Ler e verificar modelo já treinado - ok
# 2. Montar lista de códons da sequência a ser analisada - ok
# 3. Concatenar os códons da sequência analisada na ordem dos códons informados pelo modelo (listOfVarSites) - ok
# 3.1. Pegar o modelo com os códons e concatenar os códons da lista anterior - ok
# 4. Buscar no modelo se essa lista já se encontra presente
# 5. Se sim, qual cluster AGUA e a annotation associada
# 6. Se não, provável sequência nova, add ao dataset de treinamento, retreinar modelo
# 7. Genoma Detective. Treinar com 50 - 100 - 500 sequências


import sys
import dill
import time
import os

from Bio import SeqIO

import numpy as np


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


# Retorna array de sequências com id e sequência

def readSequences(filepath):
    sequences = []
    for record in SeqIO.parse(filepath, "fasta"):
        sequences.append([record.id, record.seq])
    return sequences


def nucleotideToCodon(tripleNucleotides):

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

    return str(codonNumber)


def nucleotideSequenceToCodons(sequence):

    tripleNucleotidesList = [sequence[i:i + 3]
                             for i in range(0, len(sequence), 3)]

    return list(map(lambda triple: nucleotideToCodon(triple), tripleNucleotidesList))


def collectCodonsListOfVarSites(listOfVarSites, codonSequences):
    codonsListOfVarSites = []

    for varSite in listOfVarSites:
        codonsListOfVarSites.append(codonSequences[varSite])

    return codonsListOfVarSites


def searchCollectCodonsInCodeOfClass(codeOfClass, collectCodons):
    for idx, array in enumerate(codeOfClass):
        if np.array_equal(collectCodons, array):
            return idx
    return None


def print_model(model, n):
    print("name:", model.name)
    # Posição dos códons variáveis
    print("ListOfVarSites:\n", model.ListOfVarSites)
    print("NmbOfClasses:", model.NmbOfClasses)  # Número de classes
    # Códons de cada classe
    print("CodeOfClass[:", n, "]:\n", model.CodeOfClass[:n])
    # Anotação de cada classe
    print("GroundTruth[:", n, "]:\n", model.GroundTruth[:n])
    # Id de cada sequência/classe
    print("Accession[:", n, "]:\n", model.Accession[:n])
    print("repulsion:", model.repulsion)  # Parâmetro de repulsão
    print("NmbOfClusters:", model.NmbOfClusters)  # Número de clusters
    # Espécies de cada cluster
    print("TheSpeciesOfCluster:\n", model.TheSpeciesOfCluster)
    # Resolução mínima do cluster
    print("MinClusterAcc:", model.MinClusterResolution)
    print("Resolution:", model.Resolution)  # Resolução
    # Parâmetro de repulsão do CLOPE
    print("CLOPE repulsion parameter:", model.CLOPE_repulsion)
    # Matriz de distribuição
    print("DistributionMatrix:\n", model.DistributionMatrix)
    # Classe de cada cluster
    print("TheClusterOfClass:\n", model.TheClusterOfClass)


# Função retorna array com a classe e cluster da sequência
def analyseSequence(model, sequence):
    codonSequences = nucleotideSequenceToCodons(sequence[1])

    collectCodons = collectCodonsListOfVarSites(
        model.ListOfVarSites, codonSequences)

    indexCollectCodonsInCodeOfClass = searchCollectCodonsInCodeOfClass(
        model.CodeOfClass, collectCodons)

    return [model.GroundTruth[indexCollectCodonsInCodeOfClass], model.TheClusterOfClass[indexCollectCodonsInCodeOfClass]]


def backup():
    contador = [{} for _ in range(loaded_model.NmbOfClusters)]

    # Iterar sobre cada cluster
    for target in range(loaded_model.NmbOfClusters):
        # Iterar sobre cada pclass no cluster atual
        for pclass in loaded_model.TheClusterOfClass.keys():
            # Verificar se a pclass pertence ao cluster atual
            if loaded_model.TheClusterOfClass[pclass] == target:
                # Obter o valor correspondente de model.GroundTruth[pclass]
                valor = loaded_model.GroundTruth[pclass]
                # Atualizar a contagem para o valor no dicionário do cluster atual
                if valor in contador[target]:
                    contador[target][valor] += 1
                else:
                    contador[target][valor] = 1

    print(contador)


if __name__ == "__main__":
    # Verifica se os argumentos necessários foram passados
    if len(sys.argv) < 3:
        print(
            "Uso: python3 AGUA_classificação.py <arquivo_modelo> <arquivo_sequencia> <pasta_resultados>")
        sys.exit(1)

    # Caminho do arquivo do modelo
    model_filepath = sys.argv[1]

    # Caminho do arquivo da sequência
    sequence_filepath = sys.argv[2]

    # Caminho da pasta de resultados
    result_path = os.path.join(sys.argv[3], "results.csv")

    # Recarregar o arquivo
    with open(model_filepath, "rb") as file:
        loaded_model = dill.load(file)

    # Criar o arquivo de resultados
    with open(result_path, "w") as f:
        f.write(f"id_sequencia;tipo_anotação;cluster_agua;tempo_processamento\n")

    # Ler as sequências a serem analisadas
    sequences = readSequences(sequence_filepath)

    # Contar o tempo de execução da análise de cada sequência e imprimir o resultado
    fullTime = 0
    for sequence in sequences:
        startTime = time.time()

        result_analyse = analyseSequence(loaded_model, sequence)
        print(
            f"Análise: Sequência: {sequence[0]} - Tipo: {result_analyse[0]} - Cluster: {result_analyse[1]}")

        endTime = time.time()  # Fim da medição do tempo
        executionTime = endTime - startTime  # Cálculo da duração

        executionTimeStr = str(executionTime).replace('.', ',')

        with open(result_path, "a") as f:
            f.write(
                f"{sequence[0]};{result_analyse[0]};{result_analyse[1]};{executionTimeStr}\n")

        fullTime += executionTime

        startTime = endTime = 0

        # print(f"Sequência {sequence[0]} processada em: {executionTimeStr} segundos")

    print(f"Tempo total da análise: {fullTime} segundos")
