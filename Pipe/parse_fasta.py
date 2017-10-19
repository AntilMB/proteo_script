if __name__ == '__main__':
	f = open('NCBInr_20160328.fasta', 'r')

	count = 0
	out = open('fasta_parsed/' + str(count) + '.fasta', 'w')
	key = True

	for line in f:

		if key:
			 if line.startswith('>'):
			 	key = False
			 else:
			 	continue

		if count % 20000 == 0:
			out.close()
			out = open('fasta_parsed/' + str(count) + '.fasta', 'w')

		line = line.strip()
		if line.startswith('>'):
			count += 1

		print(line, file=out)

	out.close()
	f.close()



