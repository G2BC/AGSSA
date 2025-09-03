import sys

def carregar_cabecalhos_sequencias(arquivo_fasta):
    cabecalhos = []
    with open(arquivo_fasta, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                # Remove o '>' e espaços em branco
                cabecalho = line[1:].strip()
                cabecalhos.append(cabecalho)
    return cabecalhos

def atualizar_anotacoes(arquivo_anotacoes, cabecalhos_validos):
    linhas_atualizadas = []
    with open(arquivo_anotacoes, 'r') as anotacoes_file:
        anotacoes_dict = {line.split(',')[0].strip(): line for line in anotacoes_file}

    for cabecalho in cabecalhos_validos:
        if cabecalho in anotacoes_dict:
            linhas_atualizadas.append(anotacoes_dict[cabecalho])

    # Reescreve o arquivo de anotações com as linhas atualizadas
    with open(arquivo_anotacoes, 'w') as anotacoes_file:
        anotacoes_file.writelines(linhas_atualizadas)

def main():
    if len(sys.argv) != 3:
        print("Uso: python atualizar_anotacoes.py <arquivo_sequencias> <arquivo_anotacoes>")
        sys.exit(1)

    arquivo_sequencias = sys.argv[1]
    arquivo_anotacoes = sys.argv[2]

    # Carrega os cabeçalhos das sequências
    cabecalhos_validos = carregar_cabecalhos_sequencias(arquivo_sequencias)

    # Atualiza o arquivo de anotações
    atualizar_anotacoes(arquivo_anotacoes, cabecalhos_validos)

    print(f"Anotações atualizadas salvas em {arquivo_anotacoes}")

if __name__ == "__main__":
    main()

