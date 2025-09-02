import os
import subprocess
import sys

# Função para executar um comando e retornar a saída


def run_command(command):
    result = subprocess.run(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result.stdout.strip(), result.stderr.strip(), result.returncode


if __name__ == "__main__":
    # Verifica se os argumentos necessários foram passados
    if len(sys.argv) < 3:
        print(
            "Uso: python 02_alinhar_sequencias.py <arquivo_sequencias> <arquivo_referencia>")
        sys.exit(1)

    # Caminho do arquivo específico e do arquivo de referência
    arquivo_sequencias = sys.argv[1]
    arquivo_referencia = sys.argv[2]

    # Verifica se o arquivo de sequência existe
    if os.path.exists(arquivo_sequencias):
        # Define o caminho completo para o arquivo de saída .sam
        arquivo_sam = arquivo_sequencias.replace(
            "sequencias_unicas.fasta", "mysam.sam")

        # Executa o primeiro comando para gerar o arquivo .sam
        command1 = f"minimap2 -a -x asm5 --cs --sam-hit-only --secondary=no -t 10 {arquivo_referencia} {arquivo_sequencias} -o {arquivo_sam}"
        stdout1, stderr1, returncode1 = run_command(command1)
        print(f"Comando 1: {command1}")
        print(f"stdout1: {stdout1}")
        print(f"stderr1: {stderr1}")

        # Verifica se o arquivo .sam foi gerado com sucesso
        if returncode1 == 0 and os.path.exists(arquivo_sam):
            # Define o caminho completo para o arquivo de saída sequencias_alinhadas.fasta
            arquivo_alinhado = arquivo_sequencias.replace(
                "sequencias_unicas.fasta", "sequencias_alinhadas.fasta")

            # Executa o segundo comando para gerar o arquivo sequencias_alinhadas.fasta
            command2 = f"gofasta sam toMultiAlign -s {arquivo_sam} -t 10 --reference {arquivo_referencia} > {arquivo_alinhado}"
            stdout2, stderr2, returncode2 = run_command(command2)
            print(f"Comando 2: {command2}")
            print(f"stdout2: {stdout2}")
            print(f"stderr2: {stderr2}")

            # Verifica se o arquivo sequencias_alinhadas.fasta foi gerado com sucesso
            if returncode2 == 0 and os.path.exists(arquivo_alinhado):
                print("Processamento concluído com sucesso.")
            else:
                print(
                    f"Erro ao gerar o arquivo sequencias_alinhadas.fasta para {arquivo_sequencias}")
                print(f"Saída de erro: {stderr2}")
        else:
            print(f"Erro ao gerar o arquivo .sam para {arquivo_sequencias}")
            print(f"Saída de erro: {stderr1}")
    else:
        print(f"Arquivo {arquivo_sequencias} não encontrado.")
