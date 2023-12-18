#!/usr/bin/env python3

TOTAL = 0
CONTADOR_NOTAS = 0

while CONTADOR_NOTAS < 10: #script alterado para que o número de notas a ser inserido seja 10
    TOTAL += float(input("Digite a nota: ")) #transformando o tipo de input para número flutuante #pedindo para que o usuário insira as notas e somando essas na variável TOTAL
    CONTADOR_NOTAS += 1 #impedindo com que o loop seja infinito

MEDIA = (TOTAL/10)

print (MEDIA) 

