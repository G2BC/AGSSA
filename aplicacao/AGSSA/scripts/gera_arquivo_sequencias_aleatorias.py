import sys
import random
from Bio import SeqIO

def main(input_file, num_sequences):
    # Ler todas as sequências do arquivo FASTA
    sequences = list(SeqIO.parse(input_file, "fasta"))

    # Verificar se o número de sequências solicitado é maior do que o número disponível
    if num_sequences > len(sequences):
        print(f"O arquivo contém apenas {len(sequences)} sequências. Não é possível extrair {num_sequences} sequências.")
        return

    # Selecionar N sequências aleatoriamente
    selected_sequences = random.sample(sequences, num_sequences)

    # Criar o nome do arquivo de saída
    output_file = f"sequencias_aleatorias_{num_sequences}.fasta"

    # Escrever as sequências selecionadas no arquivo de saída
    SeqIO.write(selected_sequences, output_file, "fasta")

    print(f"Arquivo '{output_file}' criado com {num_sequences} sequências aleatórias.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python script.py <caminho_arquivo_fasta> <numero_sequencias>")
    else:
        input_file = sys.argv[1]
        num_sequences = int(sys.argv[2])
        main(input_file, num_sequences)

