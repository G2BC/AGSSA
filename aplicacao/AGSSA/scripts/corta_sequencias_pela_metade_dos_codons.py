import sys
from Bio import SeqIO

def carregar_sequencias_mesmo_tamanho(arquivo_fasta):
    sequencias = list(SeqIO.parse(arquivo_fasta, "fasta"))
    tamanho_sequencias = [len(seq.seq) for seq in sequencias]
    tamanho_comum = max(set(tamanho_sequencias), key=tamanho_sequencias.count)
    
    return [seq for seq in sequencias if len(seq.seq) == tamanho_comum], tamanho_comum

def cortar_sequencias_ao_meio(sequencias, tamanho_comum):
    meio_codons = (tamanho_comum // 3) // 2
    nucleotideos_corte = meio_codons * 3
    
    sequencias_cortadas = []
    for seq in sequencias:
        metade_sequencia = seq[:nucleotideos_corte]
        sequencias_cortadas.append(metade_sequencia)
    
    return sequencias_cortadas

def salvar_sequencias_cortadas(sequencias_cortadas, arquivo_saida):
    SeqIO.write(sequencias_cortadas, arquivo_saida, "fasta")
    print(f"Arquivo '{arquivo_saida}' criado com {len(sequencias_cortadas)} sequências cortadas ao meio.")

def main():
    if len(sys.argv) != 2:
        print("Uso: python cortar_sequencias.py <arquivo_fasta>")
        sys.exit(1)
    
    arquivo_fasta = sys.argv[1]
    sequencias, tamanho_comum = carregar_sequencias_mesmo_tamanho(arquivo_fasta)
    
    if not sequencias:
        print("Nenhuma sequência encontrada com o tamanho comum.")
        return
    
    sequencias_cortadas = cortar_sequencias_ao_meio(sequencias, tamanho_comum)
    arquivo_saida = "sequencias_cortadas_ao_meio.fasta"
    salvar_sequencias_cortadas(sequencias_cortadas, arquivo_saida)

if __name__ == "__main__":
    main()

