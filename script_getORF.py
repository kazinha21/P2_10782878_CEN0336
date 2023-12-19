#!/usr/bin/env python3

import sys
import re #importanto bibliotecas necessárias

def pegar_sequencias(arquivo): #função que recebe um caminho para um arquivo multifasta e retorna todas as sequências em um dicionário 
    sequencias = {} #criar dicionário das sequências para amazená-las
    try:
        with open (arquivo, "r") as fasta: #abrir arquivo
            id_atual = "" #criar string para armazenar identificadores atuais
            for linha in fasta:
                linha = linha.rstrip() #retirar espaço e quebra de linha 
                if linha.startswith(">"):
                    id_atual = linha[1:].split()[0] #pegar o nome da sequência se essa linha for um cabeçalho
                    sequencias[id_atual] = "" #iniciar dicionário
                else: 
                    sequencias[id_atual] += linha.upper() #concatenar a linha à sequência e passar os caracteres da sequência para maiúsculo
    except FileNotFoundError:
        print (f"Arquivo não encontrado")
        exit() #mostrar erro e terminar o programa caso o arquivo não seja encontrado
    return sequencias #retornar os resultados da função

def reverso_complementar(sequencia): #função que recebe uma sequência de DNA e retorna o reverso complementar
    reverseSeq = sequencia[::-1] #pegar o reverso
    reverseSeq = reverseSeq.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c") #pegar o complementar
    return reverseSeq.upper() #retornar os resultados da função com os caracteres da sequência em maiúsculo

def maior_ORF(sequencia): #função que recebe uma sequência de DNA e retorna o maior ORF
    if len(sequencia) == 0:
        print (f" A sequência está vazia") #mostrar erro
    big_frame = 0 #criar variáveis para armazenar resultados posteriores
    big_start = 0
    big_end = 0
    big_sequence = "" 
    stop_codons = ['TAA', 'TAG', 'TGA'] #armazenar stop códons

    reverseSeq = reverso_complementar(sequencia) #utilizar a função para encontrar o reverso complementar

    for frame in range (6): #encontrar o maior ORF dentre os 6 quadros de leitura
        if frame >= 3:
            dna = reverseSeq[frame-3:] #frames maiores ou igual a 3 corresponderão aos frames negativos
        else:
            dna = sequencia[frame:] #começar de um nucleotídeo posterior para cada frame
        codons = re.findall(r".{3}", dna) #códons são grupos de 3 nucleotídeos
        frameId = frame + 1 #adicionar mais 1 ao número do frame, visto que python se inicia em 0
        
        seq = "" #criar variável para armazenar resultados posteriores
        i = 0 #número do códon
        for codon in codons:
            i += 1 #atualizar o número do códon
            if codon == "ATG" and seq == "": #se começar com ATG e se já não estiver dentro de um ORF, começar a sequência
                seq += codon
            elif codon in stop_codons and len(seq) > 0: #se estiver dentro de um ORF e for códon de parada, terminar a sequência
                seq += codon
                if len(seq) > len(big_sequence): #caso essa sequência seja maior que a anterior, atualizar os resultados
                    big_sequence = seq #armazenar a sequência
                    big_frame = frameId #armazenar o frame
                    big_start = (i*3) - len(seq) + 1 + frame%3 #o começo dessa ORF na sequência é o final da sequência menos o seu tamanho e mais um
                    big_end = (i*3) + frame%3 #o final da sequência é o número do códon vezes 3
                seq = "" #recomeçar a sequência
            elif len(seq) > 0: #se estiver dentro de um ORF, adicionar o códon na sequência
                seq +=codon 
    return big_sequence, big_frame, big_start, big_end #retornar os resultados da função

def traducao(sequencia): #função que recebe uma sequência de DNA e retorna sua tradução
    translation_table = {
        'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
        'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
        'AAT':'N', 'AAC':'N',
        'GAT':'D', 'GAC':'D',
        'TGT':'C', 'TGC':'C',
        'CAA':'Q', 'CAG':'Q',
        'GAA':'E', 'GAG':'E',
        'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
        'CAT':'H', 'CAC':'H',
        'ATT':'I', 'ATC':'I', 'ATA':'I',
        'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
        'AAA':'K', 'AAG':'K',
        'ATG':'M',
        'TTT':'F', 'TTC':'F',
        'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
        'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
        'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
        'TGG':'W',
        'TAT':'Y', 'TAC':'Y',
        'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
        'TAA':'*', 'TGA':'*', 'TAG':'*'
    } #dicionário dos aminoácidos e seus códons

    protein = "" #criar variável para armazenar resultados posteriores
    for i in range(0, len(sequencia), 3): #ir de 3 em 3 nucleotídeos
        codon = sequencia[i: i+3] #pegar códon
        if codon in translation_table:
            protein += translation_table[codon] #adicionar o aminoácido no peptídeo final
        else:
            print(f"Códon não encontrado") #mostrar erro 
    return protein #retornar os resultados da função

def main(): #função principal
    try:
        arquivo = sys.argv[1] #pegar o arquivo 
    except IndexError:
        print(f"Forneça um arquivo de entrada") #mostrar erro e terminar o programa caso o arquivo não seja fornecido
        exit()

    sequencias = pegar_sequencias(arquivo) #utilizar a função para pegar as sequências
    with open ("ORF.fna", "w") as fna, open ("ORF.faa", "w") as faa: #criar os arquivos novos

        for identificador in sequencias:
            sequencia, frame, start, end = maior_ORF(sequencias[identificador]) #utilizar e armazenar os resultados da função
            peptideo = traducao(sequencia) #utilizar a função para traduzir a proteína
            fna.write(f">{identificador}_frame{frame}_{start}_{end}\n{sequencia}\n")
            faa.write(f">{identificador}_frame{frame}_{start}_{end}\n{peptideo}\n") #salvar as informações nos arquivos

if __name__ == "__main__":
    main() #caso o programa seja executado na linha de comando, executar a função main 


