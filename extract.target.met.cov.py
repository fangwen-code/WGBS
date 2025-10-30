#!/usr/bin/env python3


import os
import sys
import gzip 
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='DNA methylation analysis workflow')
parser.add_argument('-s', '--sample', help='Sample name')
parser.add_argument('-d', '--dp', type=int, default=20, help='Cutoff based on methylation sites depth')
parser.add_argument('-c', '--cov', help='Coverage file of methylation sites from bismark')
parser.add_argument('-b', '--bed', help='The bed file of target region')
parser.add_argument('-o', '--outDir', help='The output directory')
args = parser.parse_args()

dp = args.dp
sample = args.sample
covFile = args.cov
bedFile = args.bed
outDir = args.outDir

if not os.path.exists(outDir):
    os.makedirs(outDir)

regionDict = defaultdict(list)

with open(bedFile, 'r') as f_in1:
    for rawLine in f_in1.readlines():
        arr = rawLine.strip().split('\t')
        start, end = map(int, [arr[1],arr[2]])
        regionDict[arr[0]].append([start,end])

outFile = outDir+'/'+sample+'.target.bismark.filter.cov.gz'
f_out = gzip.open(outFile, 'wb')

with gzip.open(covFile, 'rb') as f_in2:
    for line in f_in2.readlines():
        cont = line.decode()
        posArr = cont.strip().split('\t')
        chrom = posArr[0]
        if not chrom.startswith('chr'):
            chrom = 'chr'+chrom
        pos = int(posArr[1])
        depth = int(posArr[4]) + int(posArr[5])
        if chrom in regionDict:
            posLst = regionDict[chrom]
            for posReg in posLst:
                if pos>=posReg[0] & pos<=posReg[1]:
                    if depth >= dp:
                        f_out.write(cont.encode())
                        break


f_out.close()
