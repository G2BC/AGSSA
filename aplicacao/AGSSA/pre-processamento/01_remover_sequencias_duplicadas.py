import os
import sys
from Bio import SeqIO


def remove_duplicates(input_file, output_file):
    # Input
    input_obj = open(input_file, 'r')

    # Output
    output_obj = open(output_file, 'w')

    # Dictionary to store unique sequences by their sequence string
    uniq_seqs = {}

    # Iterate through input sequences
    for qry in SeqIO.parse(input_obj, 'fasta'):
        seq_str = str(qry.seq)

        # Check if this sequence string has been seen before
        if seq_str not in uniq_seqs:
            # If it's a new sequence string, save it
            uniq_seqs[seq_str] = qry

    # Making unique sequences
    final_seq = list(uniq_seqs.values())

    # Write output file
    SeqIO.write(final_seq, output_obj, 'fasta')

    # Close objects
    input_obj.close()
    output_obj.close()


if __name__ == "__main__":
    # Caminho para as sequências
    if len(sys.argv) > 1:
        arquivo_sequencias = sys.argv[1]
        print(arquivo_sequencias)

        if os.path.exists(arquivo_sequencias):
            pasta_upload = os.path.dirname(arquivo_sequencias)

            path_sequencias_filtradas = os.path.join(
                pasta_upload, 'sequencias_unicas.fasta')

            # Chamar a função para remover duplicatas e salvar o arquivo final
            remove_duplicates(arquivo_sequencias, path_sequencias_filtradas)
        else:
            print(f"Arquivo {arquivo_sequencias} não encontrado.")
    else:
        print("Informe o caminho do arquivo")
