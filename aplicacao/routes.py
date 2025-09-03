from flask import render_template, request, redirect, session, url_for, jsonify
from extensions import app, celery, aplicacao_path,AGSSA_analise_path,AGSSA_treinamento_path,UPLOADS_PATH,RESULTS_PATH
from tasks import pre_processamento_sequencias, process_analyses, process_files_and_send_email
from utils.file_utils import conta_quantidade_sequencias, ler_especies
import os, datetime, shutil, uuid
from hashlib import sha256 as _hash
from dotenv import load_dotenv

load_dotenv()




@app.route('/')
def index():
    return render_template('index.html')


@app.route('/treinamento-direto', methods=['GET', 'POST'])
def treinamento_direto():
    if request.method == 'POST':
        sequencias = request.files['sequencias']
        anotacoes = request.files['anotacoes']
        email = request.form['email']
        especie = request.form['especie']

        current_dir = os.path.abspath(os.path.dirname(__file__))
        dir_sequence = str(uuid.uuid4())

        upload_path = os.path.join(current_dir, f'{UPLOADS_PATH}/{dir_sequence}')
        result_path = os.path.join(current_dir, f'{RESULTS_PATH}/{dir_sequence}')

        os.makedirs(upload_path, exist_ok=True)
        os.makedirs(result_path, exist_ok=True)

        sequencia_path = ''
        anotacoes_path = ''

        informacoes_processamento_path = os.path.join(
            upload_path, 'informacoes_processamento.txt')

        with open(informacoes_processamento_path, 'w') as f:
            f.write(
                f'Data de Início: {datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")}\n')

        # Salvar o arquivo de sequências
        if sequencias:
            sequencia_path = os.path.join(
                upload_path, 'sequencias_treinamento.fasta')
            sequencias.save(sequencia_path)

        # Salvar o arquivo de anotações
        if anotacoes:
            anotacoes_path = os.path.join(upload_path, 'anotacoes.txt')
            anotacoes.save(anotacoes_path)

        # Copiar sequência e anotação processada para a pasta de resultados
        shutil.copy(sequencia_path, result_path)
        shutil.copy(anotacoes_path, result_path)

        quantidade_sequencias_enviadas = conta_quantidade_sequencias(
            sequencia_path)

        with open(informacoes_processamento_path, 'a') as f:
            f.write(
                f'Quantidade de sequências enviadas: {quantidade_sequencias_enviadas}\n')

        # Enviar tarefa de execução do AGSSA e envio de e-mail para o Celery
        task = process_files_and_send_email.delay(
            AGSSA_treinamento_path, sequencia_path, anotacoes_path, result_path, email, especie, informacoes_processamento_path)

        session['upload_success'] = True

        return redirect(url_for('index'))

    especies = ler_especies()
    return render_template('treinamento_direto.html', especies=especies)


@app.route('/treinamento', methods=['GET', 'POST'])
def treinamento():
    if request.method == 'POST':
        sequencias = request.files['sequencias']
        anotacoes = request.files['anotacoes']
        email = request.form['email']
        especie = request.form['especie']

        current_dir = os.path.abspath(os.path.dirname(__file__))
        dir_sequence = str(uuid.uuid4())

        upload_path = os.path.join(current_dir, f'{UPLOADS_PATH}/{dir_sequence}')
        result_path = os.path.join(current_dir, f'{RESULTS_PATH}/{dir_sequence}')

        os.makedirs(upload_path, exist_ok=True)
        os.makedirs(result_path, exist_ok=True)

        sequencia_path = ''
        anotacoes_path = ''

        informacoes_processamento_path = os.path.join(
            upload_path, 'informacoes_processamento.txt')

        with open(informacoes_processamento_path, 'w') as f:
            f.write(
                f'Data de Início: {datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")}\n')

        # Salvar o arquivo de sequências
        if sequencias:
            sequencia_path = os.path.join(
                upload_path, 'sequencias_treinamento.fasta')
            sequencias.save(sequencia_path)

        # Salvar o arquivo de anotações
        if anotacoes:
            anotacoes_path = os.path.join(upload_path, 'anotacoes.txt')
            anotacoes.save(anotacoes_path)

        # Enviar tarefa de pré-processamento para o Celery
        pre_proc_task = pre_processamento_sequencias.delay(
            sequencia_path, anotacoes_path, especie)
        sequencia_preprocessadas_path = pre_proc_task.get()

        # Copiar sequência e anotação processada para a pasta de resultados
        shutil.copy(sequencia_preprocessadas_path, result_path)
        shutil.copy(anotacoes_path, result_path)

        quantidade_sequencias_enviadas = conta_quantidade_sequencias(
            sequencia_path)
        quantidade_sequencias_processadas = conta_quantidade_sequencias(
            sequencia_preprocessadas_path)

        with open(informacoes_processamento_path, 'a') as f:
            f.write(
                f'Quantidade de sequências enviadas: {quantidade_sequencias_enviadas}\n')
            f.write(
                f'Quantidade de sequências processadas: {quantidade_sequencias_processadas}\n')

        # Enviar tarefa de execução do AGSSA e envio de e-mail para o Celery
        task = process_files_and_send_email.delay(
            AGSSA_treinamento_path, sequencia_preprocessadas_path, anotacoes_path, result_path, email, especie, informacoes_processamento_path)

        session['upload_success'] = True

        return redirect(url_for('index'))

    especies = ler_especies()
    return render_template('treinamento.html', especies=especies)


@app.route('/analise', methods=['GET', 'POST'])
def analise():
    if request.method == 'POST':
        sequencias = request.files['sequencias']
        modelo = request.files['modelo']

        current_dir = os.path.abspath(os.path.dirname(__file__))
        dir_sequence = str(uuid.uuid4())

        
        upload_path = os.path.join(
            current_dir, f'{UPLOADS_PATH}/analises/{dir_sequence}')

        os.makedirs(upload_path, exist_ok=True)

        sequencia_path = ''
        modelo_path = ''

        informacoes_processamento_path = os.path.join(
            upload_path, 'informacoes_processamento.txt')

        with open(informacoes_processamento_path, 'w') as f:
            f.write(
                f'Data de Início: {datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")}\n')

        # Salvar o arquivo de sequências
        if sequencias:
            sequencia_path = os.path.join(
                upload_path, 'sequencias_treinamento.fasta')
            sequencias.save(sequencia_path)

        # Salvar o arquivo de anotações
        if modelo:
            modelo_path = os.path.join(upload_path, 'modelo.obj')
            modelo.save(modelo_path)

        quantidade_sequencias_enviadas = conta_quantidade_sequencias(
            sequencia_path)

        with open(informacoes_processamento_path, 'a') as f:
            f.write(
                f'Quantidade de sequências analisadas: {quantidade_sequencias_enviadas}\n')

        # Enviar tarefa de execução do AGSSA e envio de e-mail para o Celery
        task = process_analyses.delay(
            AGSSA_analise_path, modelo_path, sequencia_path, upload_path, informacoes_processamento_path)

        return redirect(url_for('index'))

    return render_template('analise.html')


@app.route('/clear_session', methods=['POST'])
def clear_session():
    session.pop('upload_success', None)
    return jsonify({'status': 'session cleared'})
