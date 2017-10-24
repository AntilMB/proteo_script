# -*- coding: utf-8 -*-

from __future__ import division, absolute_import, print_function
import os
from multiprocessing import Pool, freeze_support
from os.path import isfile, join, basename
from os import listdir
import os
import re
import sys
from re import search
import argparse
from collections import defaultdict


def pooled_func(arg):
	in_fasta, out_dir, fasta = arg

	class AhoNode:
		def __init__(self):
			self.goto = {}
			self.out = []
			self.fail = None


	def aho_create_forest(patterns):
		root = AhoNode()

		for path in patterns:
			node = root
			for symbol in path:
				node = node.goto.setdefault(symbol, AhoNode())
			node.out.append(path)
		return root


	def aho_create_statemachine(patterns):
		root = aho_create_forest(patterns)
		queue = []
		for node in root.goto.itervalues():
			queue.append(node)
			node.fail = root
	
		while len(queue) > 0:
			rnode = queue.pop(0)
	
			for key, unode in rnode.goto.iteritems():
				queue.append(unode)
				fnode = rnode.fail
				while fnode is not None and key not in fnode.goto:
					fnode = fnode.fail
				unode.fail = fnode.goto[key] if fnode else root
				unode.out += unode.fail.out
	
		return root
	
	
	def aho_find_all(s, root, callback, file, prot, fasta):
		node = root
	
		for i in xrange(len(s)):
			while node is not None and s[i] not in node.goto:
				node = node.fail
			if node is None:
				node = root
				continue
			node = node.goto[s[i]]
			for pattern in node.out:
				callback(i - len(pattern) + 1, pattern, file, prot, fasta)

	
	def on_occurence(pos, patterns, file, prot, fasta):
		out = "At pos %s %s found pattern: %s" % (pos, prot, patterns)
		#print(out)
		print('\t'.join([patterns, str(pos), fasta, prot]), file=file)

	
	def read_fastaK(files):
		x = open(files)
		key = None
		seq = ""
		out = dict()
		for line in x:
			line = line.strip()
			if key is None and line.startswith(">"):
				key = line[1:]
			elif not(key is None) and line.startswith(">"):
				out[key] = seq
				seq = ""
				key = line[1:]
			else:
				seq += line.upper()
		out[key] = seq
		return(out)

	
	def replace_I_to_L(fasta):
		replaced = defaultdict(list)
		for key, value in fasta.iteritems():
			fasta[key] = value.replace('I', 'L')
			replaced[value.replace('I', 'L')].append(value)

		return replaced, fasta




	fasta_loaded = read_fastaK(in_fasta)
	replaced, fasta_loaded = replace_I_to_L(fasta_loaded)

	patterns = list(set([value for _, value in fasta_loaded.iteritems()]))
	root = aho_create_statemachine(patterns)	

	with open(out_dir + '/' + basename(fasta) + 'result.txt', 'w') as f:
			print(fasta)
			fasta_ref = read_fastaK(fasta)
			for key, value in fasta_ref.iteritems():
				aho_find_all(value, root, on_occurence, f, key, basename(fasta))

	return replaced


def make_ref(file_in, file_out):	
	if not os.path.exists(file_out):
		os.makedirs(file_out)

	f = open(file_in, 'r')	

	count = 0
	out = open(file_out + '/' + str(count) + '.fasta', 'w')
	key = True

	for line in f:
		if key:
			 if line.startswith('>'):
			 	key = False
			 else:
			 	continue

		if count % 200000 == 0:
			out.close()
			out = open(file_out + '/' + str(count) + '.fasta', 'w')

		line = line.strip()
		if line.startswith('>'):
			count += 1
			print(line, file=out)
		else:
			print(line.upper().replace('I', 'L'), file=out)

	out.close()
	f.close()



def concat_res(out_dir, replaced):
	only_txt = sorted([out_dir + '/' + f for f in listdir(out_dir) if isfile(join(out_dir+'/', f)) and search('.*?\.txt$', f)])

	with open('Pipe_result.txt', 'w') as concat_file:
		for sub_file in only_txt:
			for line in open(sub_file):
				line = line.strip()

				splited = line.split('\t')
				dif_fasta_header = splited[-1].split('\x01')
				print(splited)

				for single_header in dif_fasta_header:
					for true_pep in replaced[splited[0]]:
						print('\t'.join([true_pep] + splited[1:-1])+'\t'+single_header, file=concat_file)


def search(args):
	if not os.path.exists(args.out_dir):
		os.makedirs(args.out_dir)

	only_fasta = sorted([args.fasta_base + '/' + f for f in listdir(args.fasta_base) if isfile(join(args.fasta_base+'/', f)) and search('.*?\.fasta$', f)])
	
	pool = Pool(args.nthreads)
	#pool_arg = ['$$'.join([args.in_fasta, args.out_dir, fasta]) for fasta in only_fasta]
	pool_arg = [(args.in_fasta, args.out_dir, fasta) for fasta in only_fasta]		
	replaced = pool.map(pooled_func, pool_arg)
	# говнокод из-за особенностей pool
	replaced = replaced[0]
	concat_res(args.out_dir, replaced)


def create_base(args):
	make_ref(args.in_fasta, args.fasta_base)



if __name__ == '__main__':
	freeze_support()

	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers(dest='cmd')
	subparsers.required = True

	make_base_parser = subparsers.add_parser('create_base', help='Созадть базу для поиска')
	make_base_parser.set_defaults(_action=create_base)
	make_base_parser.add_argument('--fasta', dest='in_fasta', required=True, help='Input fasta for base')
	make_base_parser.add_argument('--base_dir', dest='fasta_base', default=None, help='parsed fasta')

	search_parser = subparsers.add_parser('search', help='Поиск по паттерну')
	search_parser.set_defaults(_action=search)
	search_parser.add_argument('--nthreads', dest='nthreads', default=20, help="nthreads")
	search_parser.add_argument('--fasta', dest='in_fasta', required=True, help='Input fasta for search')
	search_parser.add_argument('--base_dir', dest='fasta_base', required=True, help='parsed fasta')
	search_parser.add_argument('--out', dest='out_dir', required=True, help='output')

	args = parser.parse_args()
	#print(args)