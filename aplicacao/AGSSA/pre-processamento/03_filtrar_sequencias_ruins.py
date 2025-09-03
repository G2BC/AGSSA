import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

# Função para verificar se uma sequência possui mais de 30 N's consecutivos


def possui_mais_de_40_ns(sequencia):
    # Define a regex para encontrar mais de 40 'N', 'n' ou '-' seguidos
    padrao = re.compile(r'[Nn-]{41,}')
    # Verifica se o padrão existe na sequência
    return bool(padrao.search(sequencia))


def read_gene_sequences(file_genome_sequences):
    multiple_genomes = []
    with open(file_genome_sequences) as handle:
        for values in SeqIO.FastaIO.SimpleFastaParser(handle):
            sequence_id, sequence = values
            multiple_genomes.append((sequence_id, sequence))
    return multiple_genomes


if __name__ == "__main__":
    # Verifica se os argumentos necessários foram passados
    if len(sys.argv) < 2:
        print("Uso: python script.py <arquivo_sequencias>")
        sys.exit(1)

    # Caminho do arquivo de sequências
    arquivo_sequencias = sys.argv[1]

    # Verifica se o arquivo de sequência existe
    if os.path.exists(arquivo_sequencias):
        pasta_upload = os.path.dirname(arquivo_sequencias)
        path_sequencias_filtradas = os.path.join(
            pasta_upload, 'sequencias_tratadas.fasta')

        # Ler as sequências do arquivo usando a Biopython
        sequencias = read_gene_sequences(arquivo_sequencias)

        # Filtrar sequências com mais de 30 N's consecutivos
        sequencias_filtradas = []
        for seq_tuple in sequencias:
            seq_id, seq_sequence = seq_tuple
            if not possui_mais_de_40_ns(seq_sequence):
                # Criar um objeto SeqRecord com ID e sequência
                seq_record = SeqRecord(
                    Seq(seq_sequence), id=seq_id, description='')
                sequencias_filtradas.append(seq_record)

        # Salvar arquivo com sequências filtradas
        with open(path_sequencias_filtradas, "w") as novo_arquivo:
            SeqIO.write(sequencias_filtradas, novo_arquivo, "fasta")

        print(
            f"Processamento concluído com sucesso. Arquivo salvo em {path_sequencias_filtradas}")
    else:
        print(f"Arquivo {arquivo_sequencias} não encontrado.")
