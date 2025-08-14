from extensions import celery, UPLOADS_PATH, RESULTS_PATH
import os, datetime, subprocess, shutil, time
from utils.email_utils import send_email
from dotenv import load_dotenv

load_dotenv()

# Função para executar um comando e retornar a saída


# def run_command(command):
#     result = subprocess.run(
#         command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#     return result.stdout.strip(), result.stderr.strip(), result.returncode


def run_command(command):
    """Executa um comando shell e retorna stdout, stderr e returncode."""
    process = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout.decode(), stderr.decode(), process.returncode


@celery.task
def pre_processamento_sequencias(sequencia_path, anotacoes_path, especie):
    aplicacao_path = os.getenv('APLICACAO_PATH')
    referencia_path = os.path.join(
        aplicacao_path, 'database/sequencia_referencia')
    try:
        # Caminhos dos scripts
        script1 = 'AGUA/pre-processamento/01_remover_sequencias_duplicadas.py'
        script2 = 'AGUA/pre-processamento/02_alinhar_sequencias.py'
        script3 = 'AGUA/pre-processamento/03_filtrar_sequencias_ruins.py'
        script4 = 'AGUA/pre-processamento/04_filtrar_anotacoes.py'

        # Caminho do arquivo de referência
        referencia = os.path.join(referencia_path, f'{especie}.fasta')

        # Executa o primeiro script
        command1 = f'python3 {script1} {sequencia_path}'
        stdout1, stderr1, returncode1 = run_command(command1)
        if returncode1 != 0:
            raise Exception(f"Erro no script 01: {stderr1}")

        # Caminho da saída do primeiro script
        output1 = sequencia_path.replace(
            "sequencias_treinamento.fasta", "sequencias_unicas.fasta")

        # Executa o segundo script
        command2 = f'python3 {script2} {output1} {referencia}'
        stdout2, stderr2, returncode2 = run_command(command2)
        if returncode2 != 0:
            raise Exception(f"Erro no script 02: {stderr2}")

        # Caminho da saída do segundo script
        output2 = output1.replace(
            "sequencias_unicas.fasta", "sequencias_alinhadas.fasta")

        # Executa o terceiro script
        command3 = f'python3 {script3} {output2}'
        stdout3, stderr3, returncode3 = run_command(command3)
        if returncode3 != 0:
            raise Exception(f"Erro no script 03: {stderr3}")

        # Caminho da saída do terceiro script
        output3 = output2.replace(
            "sequencias_alinhadas.fasta", "sequencias_tratadas.fasta")

        # Executa o quarto script
        command4 = f'python3 {script4} {output3} {anotacoes_path}'
        stdout4, stderr4, returncode4 = run_command(command4)
        if returncode4 != 0:
            raise Exception(f"Erro no script 04: {stderr4}")

        return output3
    except Exception as e:
        raise e


@celery.task
def process_files_and_send_email(agua_path, sequencia_path, anotacoes_path, result_path, email, especie, informacoes_path):
    try:
        command = f'python3 {agua_path} {sequencia_path} {anotacoes_path} {result_path} {especie}'
        os.system(command)

        with open(informacoes_path, 'a') as f:
            f.write(
                f'Data de Término: {datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")}\n')

        # Compactar os arquivos no diretório result_path
        zip_filename = 'resultados.zip'
        zip_filepath = os.path.join(result_path, zip_filename)

        send_email(email, zip_filepath, informacoes_path)
        return "Done"
    except Exception as e:
        raise e


@celery.task
def process_analyses(agua_path, modelo_path, sequencias_path, result_path, informacoes_path):
    try:
        command = f'python3 {agua_path} {modelo_path} {sequencias_path} {result_path}'
        os.system(command)

        with open(informacoes_path, 'a') as f:
            f.write(
                f'Data de Término: {datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")}\n')

        return "Done"
    except Exception as e:
        raise e



def remove_directory(directory_path):
    try:
        now = time.time()
        one_week_ago = 7 * 24 * 60 * 60  
        

        if os.path.isdir(directory_path):
            info = os.stat(directory_path)
            age_seconds = now - info.st_mtime
        else:
            raise ValueError("The path provided is not a directory")

        if age_seconds > one_week_ago:
            shutil.rmtree(directory_path)
        

    except Exception as e:
        raise ValueError(f"Error removing directory: {e}")
    
@celery.task
def remove_old_directorys():
    for folder in os.listdir(UPLOADS_PATH):
        folder_fullpath = os.path.join(UPLOADS_PATH, folder)
        remove_directory(folder_fullpath)

    for folder in os.listdir(RESULTS_PATH):
        folder_fullpath = os.path.join(RESULTS_PATH, folder)
        remove_directory(folder_fullpath)