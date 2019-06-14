#!/usr/bin/env python

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

import sys
import os
from copy import deepcopy
import re

from fasta import Fasta

class InfernalRecord:
    '''
    class for record of infernal output.
    version: 1.1.2
    cmd: cmscan -Z 1028.452 --cut_ga --rfam --nohmmonly --tblout genome.tblout \
         -o genome.cmscan --fmt 2 --clanin Rfam.clanin Rfam.cm $genome
    format: --fmt 2, the coordinate is 1-based
    '''
    def __init__(self, x):
        '''
        @parameter x: a line of infernal output
        '''
        x = re.split(r'\s{2,}', x)   # the f**king seperation
        tmp_list = []
        for i in x[:-1]:
            if ' ' in i:
                tmp_list.extend(i.split())
            else:
                tmp_list.append(i)
        tmp_list.append(x[-1])
        self.target = tmp_list[1]
        self.qname = tmp_list[3]  # query name
        if tmp_list[11] == '+':
            self.qstart = int(tmp_list[9]) - 1  # query start, 1-based to 0-based
            self.qend = int(tmp_list[10])
        else:
            self.qstart = int(tmp_list[10]) - 1  # query start, 1-based to 0-based
            self.qend = int(tmp_list[9])
        self.qstrand = tmp_list[11]
        self.trunc = tmp_list[12]
        self.score = float(tmp_list[16])
        self.E = float(tmp_list[17])
        self.inc = tmp_list[18]
        self.olp = tmp_list[19]
        self.description = tmp_list[-1].strip('- ')


class Infernal:
    '''
    class for infernal: http://eddylab.org/infernal/
    version: 1.1.2
    cmd: cmscan -Z 1028.452 --cut_ga --rfam --nohmmonly --tblout genome.tblout \
         -o genome.cmscan --fmt 2 --clanin Rfam.clanin Rfam.cm $genome
    format: --fmt 2
    '''
    def __init__(self, tblout):
        '''
        @parameter tblout: tblout of infernal
        '''
        self.filename = os.path.abspath(os.path.realpath(tblout))
        self.records = list(self.__to_record())

    def __to_record(self):
        '''
        transform the output of infernal to interval format as InfernalRecord
        0-based coordinate
        '''
        try:
            fin = open(self.filename, 'r')
        except IOError as err:
            logger.error(str(err))
        else:
            for line in fin:
                line = line.strip()
                if line.startswith('#'):
                    continue
                yield InfernalRecord(line)

    def filter(self, deoverlap=True, remove_cmtrunc=True):
        '''
        '''
        if deoverlap:
            self.records = [record for record in self.records if record.olp != '=']
        if remove_cmtrunc:
            self.records = [record for record in self.records if record.trunc == 'no']

    def to_bed(self, output_prefix):
        '''
        transform aragorn_out to bed, which contains the following filed:
        rname, start, end, name, score, strand
        '''
        try:
            count = 1
            with open('%s.bed' %output_prefix, 'w') as fout:
                for record in self.records:
                    rname = record.qname
                    start = str(record.qstart)
                    end = str(record.qend)
                    name = 'infernal%s' %count + '-' + record.target
                    score = str(record.score)
                    strand = record.qstrand
                    fout.write('\t'.join([rname, start, end, name, score, strand]) + '\n')
                    count += 1
        except IOError as err:
            logger.error('Error in transforming aragorn_out to bed: %s' %str(err))


class AragornRecord:
    '''
    class for record of aragorn output.
    version: 1.2.38
    cmd: aragorn -v -d -l -t -gcstd -seq -br -fasta -o sf.aragorn.out $SF_genome
    '''
    def __init__(self, rname, start, end, strand, name, seq):
        '''
        @parameter rname: contig name
        @parameter start: start coordinate (0-based)
        @parameter end: end
        @parameter strand: strand
        @parameter name: name of tRNA
        @parameter seq: sequence of tRNA (a Sequence instance)
        '''
        self.rname = rname
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name
        self.seq = seq


class Aragorn:
    '''
    class for aragorn: http://mbio-serv2.mbioekol.lu.se/ARAGORN/
    version: 1.2.38
    cmd: aragorn -v -d -l -t -gcstd -seq -br -fasta -o sf.aragorn.out $SF_genome
    '''
    def __init__(self, aragorn_out, genome):
        '''
        @parameter aragorn_out: output of aragorn
        @genome: genome sequences that used to predict (will be used to get the strand info)
        '''
        self.filename = os.path.abspath(os.path.realpath(aragorn_out))
        self.genome = Fasta(genome)
        self.records = self.__to_record()

    def __to_record(self):
        '''
        transform the output of aragorn to interval format as AragornRecord
        0-based coordinate
        '''
        try:
            fin = open(self.filename)
        except IOError as err:
            logger.error(str(err))
        else:
            for line in fin:
                line = line.strip()

                if line in self.genome.keys:
                    rname = line

                if line.startswith('>'):
                    line = line.split()
                    seq_id = line[0].replace('>', '')
                    if 'c' in line[1]:
                        rc = True
                    else:
                        rc = False
                    start, end = line[1].strip('c[]').split(',')
                    start = int(start)-1; end = int(end)
                    strand = '-' if rc else '+'
                    seq = self.genome.get_seq(rname, start, end, rc)
                    yield AragornRecord(rname, start, end, strand, seq_id, seq)

    def to_bed(self, output_prefix):
        '''
        transform aragorn_out to bed6, which contains the following filed:
        rname, start, end, name, seq, strand
        '''
        try:
            with open('%s.bed' %output_prefix, 'w') as fout:
                for record in self.records:
                    rname = record.rname
                    start = str(record.start)
                    end = str(record.end)
                    strand = record.strand
                    name = record.name
                    fout.write('\t'.join([rname, start, end, name, str(0), strand]) + '\n')
        except IOError as err:
            logger.error('Error in transforming aragorn_out to bed: %s' %str(err))


class tRNAscanSERecord:
    '''
    class for record of tRNAscan-SE output.
    version: 2.0
    cmd: tRNAscan-SE --thread $PPN -E -o sf.tRNAscan.out -f sf.tRNAscan.structures -m sf.tRNAscan.stats -b sf.tRNAscan.bed -a sf.tRNAscan.fa -l sf.tRNAscan.log --detail -d -y $SF_genome
    '''
    def __init__(self, x):
        '''
        @parameter x: a line of tRNAscan-SE output
        '''
        x = x.split('\t')
        self.rname = x[0]
        tmp1, tmp2 = int(x[2]), int(x[3])
        if tmp1 < tmp2:
            self.strand = '+'
            self.start = tmp1 - 1
            self.end = tmp2
        else:
            self.strand = '-'
            self.start = tmp2 - 1
            self.end = tmp1
        self.name = 'tRNA-' + x[4] + x[5]
        self.inf_score = float(x[8])
        self.pseudo = True if 'pseudo' in x[-1] else False
        if int(x[6]) == 0:
            self.blockCount = 1
            self.blockSizes = (self.end-self.start,)
            self.blockStarts = (0,)
        else:
            self.blockCount = 2
            tmp3, tmp4 = int(x[6]), int(x[7])
            self.blockSizes = (min(tmp3, tmp4)-min(tmp1, tmp2), max(tmp1, tmp2)-max(tmp3, tmp4),)
            self.blockStarts = (0,max(tmp1, tmp2)-max(tmp3, tmp4)+1,)


class tRNAscanSE:
    '''
    a class for tRNAscan-SE
    version: 2.0
    cmd: tRNAscan-SE --thread $PPN -E -o sf.tRNAscan.out -f sf.tRNAscan.structures -m sf.tRNAscan.stats -b sf.tRNAscan.bed -a sf.tRNAscan.fa -l sf.tRNAscan.log --detail -d -y $SF_genome
    '''
    def __init__(self, tRNAscan_out):
        '''
        @parameter tRNAscan_out: output of tRNAscan-SE (the file of -o)
        '''
        self.filename = os.path.abspath(os.path.realpath(tRNAscan_out))
        self.records = list(self.__to_record())

    def __to_record(self):
        '''
        transform the output of aragorn to interval format as AragornRecord
        0-based coordinate
        '''
        try:
            fin = open(self.filename)
        except IOError as err:
            logger.error(str(err))
        else:
            for line in fin:
                line = line.strip()
                if line.startswith('Sequence'):
                    continue
                if line.startswith('Name'):
                    continue
                if line.startswith('-'):
                    continue
                yield tRNAscanSERecord(line)

    def filter(self, inf_score=0, remove_pseudo=True, remove_intron=True):
        '''
        '''
        self.records = [record for record in self.records if record.inf_score > inf_score]
        if remove_pseudo:
            self.records = [record for record in self.records if record.pseudo == False]
        if remove_intron:
            self.records = [record for record in self.records if record.blockCount == 1]

    def to_bed(self, output_prefix):
        '''
        transform tRNAscan-SE to bed12, which contains the following filed:
        rname, start, end, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
        '''
        try:
            with open('%s.bed' %output_prefix, 'w') as fout:
                for record in self.records:
                    rname = record.rname
                    start = str(record.start)
                    end = str(record.end)
                    name = record.name
                    score = str(record.inf_score)
                    strand = record.strand
                    thickStart = start
                    thickEnd = start
                    itemRgb = '0'
                    blockCount = str(record.blockCount)
                    blockSizes = ','.join(map(str, record.blockSizes)) + ','
                    blockStarts = ','.join(map(str, record.blockStarts)) + ','

                    fout.write('\t'.join([rname, start, end, name, score, strand, thickStart, thickEnd,
                                          itemRgb, blockCount, blockSizes, blockStarts]) + '\n')
        except IOError as err:
            logger.error('Error in transforming tRNAscan-SE to bed: %s' %str(err))


def test():
    aragorn_out = '/home/niuyw/Project/tRNA_promoter_190524/aragorn/sf.aragorn.out'
    genome = '/home/niuyw/RefData/Spodoptera_frugiperda/GCA_002213285.1_ASM221328v1_genomic.fa'
    a = Aragorn(aragorn_out, genome)
    #a.to_bed('sf.aragorn')

    # infernal
    infernal_out = '/home/niuyw/Project/tRNA_promoter_190524/infernal/sf.genome.tblout'
    b = Infernal(infernal_out)
    b.filter()
    b.to_bed('/home/niuyw/Project/tRNA_promoter_190524/infernal/sf.genome.filter')

    # tRNAscan
    #tRNAscan_out = '/home/niuyw/Project/tRNA_promoter_190524/tRNAscan-SE/sf.tRNAscan.out'
    #c = tRNAscanSE(tRNAscan_out)
    #c.filter()
    #c.to_bed('/home/niuyw/Project/tRNA_promoter_190524/tRNAscan-SE/sf.tRNAscan.filter')

if __name__ == '__main__':

    try:
        test()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(1)

