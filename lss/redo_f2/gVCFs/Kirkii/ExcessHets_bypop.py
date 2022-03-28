#!/usr/bin/python
#usage: python .py <vcf_file> <pop_file>


import sys


def filter_all_hets(list_of_GTs):
    GTs = [s[0:3] for s in list_of_GTs]
#    print(len(set(GTs)))
#    print(*GTs, sep = "\t")
    if len(set(GTs)) == 1:
         if GTs[0] == str("0/1"):
            return True
         else: return False
    else: return False


##read line pop file with PHRED threshhold
pops = []
indices = []
with open(sys.argv[2], "r") as handle:
    for line in handle:
        line = line.rstrip('\n').split('\t')
        pops.append(line)

 
##read in vcf file
filter_bool = True

with open(sys.argv[1], "r") as handle:
    for line in handle:
        if line[0] == "#": #header
            temp_line = line[0:8]
#            print(temp_line)
            if temp_line == "##FILTER":
                if filter_bool:
                    filter_bool = False
                    print('##FILTER=<ID=ExcessHetsByPop,Description="ExcessHets_from_Conover_et_al.">')
                print(line.strip())
            elif temp_line[0:6] == "#CHROM":
#                print('##FILTER=<ID=ExcessHetsByPop,Description="ExcessHets_from_Conover_et_al.">')
                print(line.strip())
                #Do stuff with pop_file to get indices of indiv in each pop
                line = line.strip().split('\t')
                for i in pops:
                    temp_index = []
                    for j in i:
                        temp_index.append(line.index(j))
                    indices.append(temp_index)
#                    print(*temp_index, sep = '\t')
                #print(*line, sep = '\t') 
            else: print(line.strip())
        else: 
            line = line.strip().split('\t')
            change = []
            for ind in indices:
                pop = [line[j] for j in ind]
                if filter_all_hets(pop):
                    change += ind
            if change:
                formats = line[8].split(':')
                if 'FT' not in formats:
                    formats.append("FT")
                    line[8] += ":FT"
                    for i, n in enumerate(line):
                        if i > 8:
                            line[i] += ":."
                f = formats.index('FT')
                for i, n in enumerate(line):
                    if i in change:
                        n = n.split(':')
                        if len(n[f]) == 1:
                            n[f] = n[f].replace('.', '') + "ExcessHetsByPop"
                        else: n[f] = n[f] + ";ExcessHetsByPop"
                        n = ':'.join(n)
                        line[i] = n
            print(*line, sep = '\t')


