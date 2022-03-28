#!/usr/bin/python

import sys

infile = sys.argv[1]

def count_hets(snp_line, queries):
    total_count = 0
    hets_count = 0
    for i in queries:
        if snp_line[i] != "./.":
            total_count += 1
        if snp_line[i] == "0/1":
            hets_count += 1
    return(hets_count,total_count)

AD1W = [9] + list(range(20, 29))
AD1L1 = list(range(29, 39))
AD1L2 = list(range(39, 49))
AD1D = list(range(49, 59))
AD2L1 = list(range(59, 69))
AD2L2 = list(range(69, 79))
AD2D = list(range(79, 89))
AD3 = list(range(89, 96))
AD4 = list(range(96, 102))
AD5 = list(range(102, 107))
AD6 = list(range(107, 113))
AD7 = list(range(113, 115))
A1W = list(range(119, 124))
A1D = list(range(124, 134))
A2D = list(range(10, 20))
D5 = list(range(115, 119))

pops = [AD1W, AD1L1, AD1L2, AD1D, AD2L1, AD2L2, AD2D, AD3, AD4, AD5, AD6, AD7, A1W, A1D, A2D, D5]

first_line = ["CHROM", "POS", "AD1W", "AD1L1", "AD1L2", "AD1D", "AD2L1", "AD2L2", "AD2D", "AD3", "AD4", "AD5", "AD6", "AD7", "A1W", "A1D", "A2D", "D5"]
print(*first_line, sep="\t")

with open(infile) as handle:
    for line in handle:
        if line[0] != "#":
            line = line.strip().split("\t")
#            for i in pops:
#                print(str(i))
#                for j in i:
#                    print(line[j])

            if len(line[3]) == 1 & len(line[4]) == 1:
                output = line[0:2]
                for h in range(9, len(line)):
                    line[h] = line[h][:3]
                for i in pops:
                    hets, total = count_hets(line, i)
                    output.append(hets)
                    output.append(total)
                print(*output, sep = "\t")
