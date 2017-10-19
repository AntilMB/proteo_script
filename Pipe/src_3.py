# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function, unicode_literals
import os
from multiprocessing import Pool
from os.path import isfile, join, basename
from os import listdir
import re
import sys
from re import search
import argparse


class AhoNode:
	''' Вспомогательный класс для построения дерева
	'''
	def __init__(self):
		self.goto = {}
		self.out = []
		self.fail = None


def aho_create_forest(patterns):
	'''Создать бор - дерево паттернов
	'''
	root = AhoNode()

	for path in patterns:
		node = root
		for symbol in path:
			node = node.goto.setdefault(symbol, AhoNode())
		node.out.append(path)
	return root


def aho_create_statemachine(patterns):
	'''Создать автомат Ахо-Корасика.
	Фактически создает бор и инициализирует fail-функции
	всех узлов, обходя дерево в ширину.
	'''
	# Создаем бор, инициализируем
	# непосредственных потомков корневого узла
	root = aho_create_forest(patterns)
	queue = []
	for node in root.goto.items():
		queue.append(node)
		node.fail = root

	# Инициализируем остальные узлы:
	# 1. Берем очередной узел (важно, что проход в ширину)
	# 2. Находим самую длинную суффиксную ссылку для этой вершины - это и будет fail-функция
	# 3. Если таковой не нашлось - устанавливаем fail-функцию в корневой узел
	while len(queue) > 0:
		rnode = queue.pop(0)

		for key, unode in rnode.goto.items():
			queue.append(unode)
			fnode = rnode.fail
			while fnode is not None and key not in fnode.goto:
				fnode = fnode.fail
			unode.fail = fnode.goto[key] if fnode else root
			unode.out += unode.fail.out

	return root


def aho_find_all(s, root, callback, file, prot):
	'''Находит все возможные подстроки из набора паттернов в строке.
	'''
	node = root

	for i in xrange(len(s)):
		while node is not None and s[i] not in node.goto:
			node = node.fail
		if node is None:
			node = root
			continue
		node = node.goto[s[i]]
		for pattern in node.out:
			callback(i - len(pattern) + 1, pattern, file, prot)


############################
# Демонстрация работы алгоритма
def on_occurence(pos, patterns, file, prot):
	out = "At pos %s %s found pattern: %s" % (pos, prot, patterns)
	print(out)


def read_fasta(files):
    x = open(files)
    key = None
    seq = ""
    out = dict()
    for line in x:
        if key is None and line.startswith(">"):
            key = line.strip()[1:]
        elif not(key is None) and line.startswith(">"):
            out[key] = seq
            seq = ""
            key = line.strip()[1:]
        else:
            seq += line.strip()
    out[key] = seq
    return(out)

def read_fastaK(files):
    x = open(files)	
    key = -1
    seq = ""
    out = dict()
    for line in x:
        if key == -1 and line.startswith(">"):
            key = 0
        elif not(key == -1) and line.startswith(">"):
            out[key] = seq
            seq = ""
            key += 1
        else:
            seq += line.strip()
    out[key] = seq
    return(out)


if __name__ == '__main__':

	print('read fasta pattern')
	patterns = [value for _, value in read_fastaK('GSSP.fasta').items()]
	root = aho_create_statemachine(patterns)

	only_fasta = sorted(['fasta_parsed/' + f for f in listdir('fasta_parsed') if isfile(join('fasta_parsed/', f)) and search('.*?\.fasta$', f)])


	with open('result.txt', 'w') as f:
		for file in only_fasta[1:4]:
			print(file)
			fasta_ref = read_fastaK('fasta_parsed/0.fasta')

			for key, value in fasta_ref.items():
				aho_find_all(value, root, on_occurence, file, key)
