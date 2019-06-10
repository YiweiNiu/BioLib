#!/usr/bin/env python

'''
for python3
'''

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

import sys
import os
from os.path import getmtime, getsize
import argparse
from copy import deepcopy
from collections import OrderedDict
from itertools import groupby

# self module
from stats import mean, median


class Sequence:
    '''
    class for sequences
    '''
    __slots__ = ['name', 'seq', 'long_name', 'len', 'type']
    def __init__(self, name=None, seq=None, long_name=None, seq_type='DNA'):
        '''
        @parameter name: name of sequence
        @parameter seq: sequence
        @parameter seq_type: type of sequence, one of DNA/RNA/Protein
        '''
        self.name = name
        self.seq = seq
        self.long_name = long_name if long_name else name
        self.len = len(seq)
        self.type = seq_type

    def __str__(self):
        return '>' + self.long_name + '\n' + '\n'.join(self.seq[i:i+60] for i in range(0, len(self.seq), 60))

    def __repr__(self):
        return 'Sequence("%s")' % (self.name)

    def __eq__(self, x):
        '''self == x
        '''
        if not isinstance(x, type(self)):
            raise TypeError('Object provided is not a Sequence instance.')
        return self.seq.lower() == x.seq.lower()

    def __ne__(self, x):
        '''self != x
        '''
        if not isinstance(x, type(self)):
            raise TypeError('Object provided is not a Sequence instance.')
        return self.seq.lower() != x.seq.lower()

    def __lt__(self, x):
        '''self < x
        '''
        if not isinstance(x, type(self)):
            raise TypeError('Object provided is not a Sequence instance.')
        return self.len < x.len

    def __le__(self, x):
        '''self <= x
        '''
        if not isinstance(x, type(self)):
            raise TypeError('Object provided is not a Sequence instance.')
        return self.len <= x.len

    def __gt__(self, x):
        '''self > x
        '''
        if not isinstance(x, type(self)):
            raise TypeError('Object provided is not a Sequence instance.')
        return self.len > x.len

    def __ge__(self, x):
        '''self >= x
        '''
        if not isinstance(x, type(self)):
            raise TypeError('Object provided is not a Sequence instance')
        return self.len >= x.len

    def __add__(self, x):
        '''self + x
        '''
        if isinstance(x, type(self)):
            if self.type != x.type:
                raise Exception('Add operation only same type sequences')
            return self.__class__(self.name+'_add_'+x.name, self.seq+x.seq)
        elif isinstance(x, str):
            self.seq += x
            return self.__class__(self.name+'_add_'+x.name, self.seq+x)
        else:
            raise Exception('Add operation only supports Sequence object or string')

    def GC(self):
        '''Return nubmer of G+C of seq as a int
        '''
        if self.type == 'Protein':
            raise Exception('Protein sequence does not support count G+C')
        g = self.seq.count('G')
        g += self.seq.count('g')
        c = self.seq.count('C')
        c += self.seq.count('c')
        return g+c

    def AT(self):
        '''Return number of A+T(U) of seq as a int
        '''
        a = self.seq.count('A')
        a += self.seq.count('a')
        if self.type == 'DNA':
            t = self.seq.count('T')
            t += self.seq.count('t')
            return a+t
        elif self.seq == 'RNA':
            u = self.seq.count('U')
            u += self.seq.count('u')
            return a+u
        else:
            raise Exception('Protein sequence does not support count A+T(U)')

    def gap(self):
        '''Return nubmer of N of seq as a int
        '''
        n = self.seq.count('N')
        n += self.seq.count('n')
        return n

    def sub(self, start=None, end=None):
        '''
        1-based corrdinate
        '''
        return self.__class__(self.name + '[%s, %s]' %(start, end), self.seq[start-1:end], 
                              self.long_name, self.type)

    def lower(self):
        '''lower case
        '''
        return self.__class__(self.name, self.seq.lower(), self.long_name, self.type)

    def upper(self):
        '''upper case
        '''
        return self.__class__(self.name, self.seq.upper(), self.long_name, self.type)

    def reverse(self):
        '''
        reverse the seq
        '''
        return self.__class__(self.name + '/r', self.seq[::-1], self.long_name, self.type)

    def revcomp(self):
        '''
        reverse complement the seq (DNA/RNA), return 5' to 3'
        '''
        seq_rv = self.seq[::-1]
        if self.type == 'DNA':
            seq_rc = seq_rv.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))
            return self.__class__(self.name + '/rc', seq_rc, self.long_name, self.type)
        elif self.type == 'RNA':
            seq_rc = seq_rv.seq.translate(str.maketrans('ACGUacguRYMKrymkVBHDvbhd', 'UGCAugcaYRKMyrkmBVDHbvdh'))
            self.name = self.name + '/rc'
            return self.__class__(self.name + '/rc', seq_rc, self.long_name, self.type)
        else:
            raise Exception('Method revcomp only supports DNA and RNA sequences.')

    def __neg__(self):
        '''reverse complement
        '''
        return self.revcomp()

    def match(self, x):
        '''
        short sequence compare
        @parameter x: a seq instance

        @return: match 1, reverse match 2, not match 0.
        '''
        if not isinstance(x, type(self)):
            raise TypeError('Object provided is not a Sequence instance.')
        seq1 = self.upper()
        seq1_rc = self.revcomp().upper()
        seq2 = x.upper()

        if seq1 == seq2:
            return '+'
        elif seq1_rc == seq2:
            return '-'
        else:
            return None


class Faidx:
    '''
    A python implementation of samtools faidx FASTA indexing
    faidx: http://www.htslib.org/doc/faidx.html
    '''
    def __init__(self, filename):
        '''
        init a faidx class

        @parameter finename: file name of the fasta file
        '''
        self.filename = os.path.abspath(os.path.realpath(filename))

        if filename.lower().endswith('.gz'):
            import gzip
            self._file_opener = gzip.open
        else:
            self._file_opener = open

        self.idxname = filename + '.fai'
        self.idx = OrderedDict()

        if os.path.exists(self.idxname) and \
           getsize(self.idxname) > 0 and \
           getmtime(self.idxname) >= getmtime(self.filename):
            self.read_idx()
        elif os.path.exists(self.idxname) and \
             getsize(self.idxname) > 0 and \
             getmtime(self.idxname) < getmtime(self.filename):
            logger.warning('Index file {0} is older than FASTA file {1}'.format(self.idxname, self.filename))
            self.read_idx()
        else:
            self.build_idx()
            self.read_idx()

    def read_idx(self):
        '''
        read index file
        '''
        try:
            with open(self.idxname, 'r') as fin:
                for line in fin:
                    line = line.strip().split('\t')
                    # seq_id: seq_len, offset, linebases, linewidth
                    self.idx[line[0]] = list(map(int, line[1:]))
        except:
            logger.error('Error in reading fasta index file.')

    def build_idx(self):
        '''
        build index file
        '''
        try:
            with self._file_opener(self.filename, 'rb') as fin:  # must be 'rb' mode
                with open(self.idxname, 'w') as fout:
                    tmp_dict = OrderedDict()
                    for line in fin:
                        line = line.lstrip()
                        if line.startswith(b'>'):
                            if 'seq_id' in vars():
                                if len(tmp_len) > 2:
                                    raise Exception('Line length of contig %s is not consistent!' %seq_id)
                                tmp_dict[seq_id][2] = max(tmp_len)
                                tmp_dict[seq_id][3] = max(tmp_width)
                            seq_id = line.split()[0].replace(b'>', b'')
                            tmp_dict[seq_id] = [0, fin.tell(), 0, 0]
                            tmp_len = set()
                            tmp_width = set()
                        else:
                            x = len(str(line.rstrip(), 'utf-8'))
                            tmp_dict[seq_id][0] += x
                            tmp_len.add(x)
                            tmp_width.add(len(line))
                    # the last seq
                    if len(tmp_len) > 2:
                        raise Exception('Line length of contig %s is not consistent!' %seq_id)
                    tmp_dict[seq_id][2] = max(tmp_len)
                    tmp_dict[seq_id][3] = max(tmp_width)

                    for seq_id in tmp_dict:
                        fout.write('\t'.join(list(map(str, [str(seq_id, 'utf-8')]+tmp_dict[seq_id]))) + '\n')

        except OSError as err:
            logger.error(str(err))

    def fetch(self, name, start, end):
        '''
        fetch seq from fasta using interval [start, end]
        1-based region, end-inclusive
        '''
        assert start == int(start)
        assert end == int(end)

        try:
            i = self.idx[name]
        except KeyError:
            raise KeyError('Requested rname {0} does not exist!'.format(name))

        if start > end:
            raise Exception('Requested end coordinate must be greater than start coordinate.')

        if start < 1:
            raise Exception('Requested start coordinate must be greater than 0.')

        if end > self.idx[name][0]:
            raise Exception('Requested end coordinate outside of {0}.'.format(name))

        size, offset, linebases, linebytes = self.idx[name]
        linediff = linebytes - linebases

        # get offset to fetch
        start -= 1
        offset += int(start / linebases) * linebytes + start % linebases
        realsize = end-start  # 1-based, inclusive
        if end > linebases:
            # if over a line, then add a byte for '\n'
            if linebases in range(start % linebases, end % linebases + linebases):
                bytesize = int(realsize / linebases) * linebytes + realsize % linebases + linediff
            else:
                bytesize = int(realsize / linebases) * linebytes + realsize % linebases
        else:
            bytesize = int(realsize / linebases) * linebytes + realsize % linebases

        try:
            fin = self._file_opener(self.filename, 'rb')
        except IOError as err:
            logger.error(str(err))
        else:
            fin.seek(offset, 0)
            seq = fin.read(bytesize).replace(b'\n', b'')
            seq = str(seq, 'utf-8')
            #print (seq)
            name = '%s:%s-%s' %(name, start+1, end)
            return Sequence(name, seq)


class Fasta:
    '''
    '''
    def __init__(self, filename=None, seq_type='DNA'):
        self.filename = os.path.abspath(os.path.realpath(filename))
        self.type = seq_type
        if filename.lower().endswith('.gz'):
            import gzip
            self._file_opener = gzip.open
        else:
            self._file_opener = open
        self.faidx = Faidx(self.filename)
        self.records = dict([(record.name, record) for record in self.__to_record()])
        self.len = len(self.faidx.idx)


    def __to_record(self):
        '''
        convert fasta file to record
        modified from: https://www.biostars.org/p/710/#120760
        '''
        try:
            fin = self._file_opener(self.filename, 'rt')
        except IOError as err:
            logger.error(str(err))
        else:
            faiter = (x[1] for x in groupby(fin, lambda line: line[0] == ">"))
            for header in faiter:
                headerStr = header.__next__()
                long_name = headerStr.strip().replace('>', '')
                name = long_name.split()[0]
                seq = "".join(s.strip() for s in faiter.__next__())
                yield Sequence(name, seq, long_name, self.type)

    def __repr__(self):
        return 'Fasta("%s")' % (self.filename)

    def __len__(self):
        return self.len

    def __iter__(self):
        for rname in self.records:
            yield self.records[rname]

    def __report(self, record, fh, reverse=False, revcomp=False, lower=False, upper=False, 
                 name_only=False, seq_only=False):
        if reverse:
            record = record.reverse()
        elif revcomp:
            record = record.revcomp()
        if lower:
            record = record.lower()
        elif upper:
            record = record.upper()
        if name_only:
            fh.write(record.name + '\n')
        elif seq_only:
            fh.write(record.seq + '\n')
        else:
            fh.write(str(record) + '\n')

    def report(self, outfilename='-', minlen=0, maxlen=-1, reverse=False, revcomp=False,
               lower=False, upper=False, name_only=False, seq_only=False):
        '''reporting fasta
        '''
        try:
            if outfilename == '-':
                fout = sys.stdout
            elif outfilename.lower().endswith('.gz'):
                import gzip
                fout = gzip.open(outfilename, 'wt')
            else:
                fout = open(outfilename, 'wt')
        except IOError as err:
            logger.error(str(err))
        else:
            for rname in self.records:
                record = self.records[rname]
                if record.len <= minlen:
                    continue
                if maxlen == -1:
                    self.__report(record, fout, reverse, revcomp, lower, upper, name_only, seq_only)
                else:
                    if record.len >= maxlen:
                        continue
                    self.__report(record, fout, reverse, revcomp, lower, upper, name_only, seq_only)
            fout.close()

    def subseq(self, name, start, end, rc=False):
        '''
        sub sequence from fasta using interval [start, end]
        1-based region, end-inclusive
        If rc is set, reverse complement will be returned.
        '''
        seq = self.faidx.fetch(name, start, end)
        if rc:
            return seq.revcomp()
        else:
            return seq

    def GC(self):
        '''%gc content of fasta, as a float
        '''
        sum_len = sum(self.faidx.idx[x][0] for x in self.faidx.idx)
        gc = sum(self.records[rname].GC() for rname in self.records)
        return gc / sum_len

    def AT(self):
        '''%at(u) content of fasta, as a float
        '''
        sum_len = sum(self.faidx.idx[x][0] for x in self.faidx.idx)
        at = sum(self.records[rname].AT() for rname in self.records)
        return at / sum_len

    def gap(self):
        '''%gap content of fasta, as a float
        '''
        sum_len = sum(self.faidx.idx[x][0] for x in self.faidx.idx)
        gap = sum(self.records[rname].gap() for rname in self.records)
        return gap / sum_len

    def N50(self, n=0.5):
        '''
        '''
        tmp_len = [self.faidx.idx[x][0] for x in self.faidx.idx]
        tmp_len = sorted(tmp_len, reverse=True)
        N50_pos = sum(tmp_len) * n

        tmp_sum = 0
        for value in tmp_len:
            tmp_sum += value
            if N50_pos <= tmp_sum:
                N50 = value
                break
        return N50

    def stats(self, outfilename='-', each=False):
        '''
        simple stat of fasta: 
        '''
        try:
            if outfilename == '-':
                fout = sys.stdout
            elif outfilename.lower().endswith('.gz'):
                import gzip
                fout = gzip.open(outfilename, 'wt')
            else:
                fout = open(outfilename, 'wt')
        except IOError as err:
            logger.error(str(err))
        else:
            tmp_len = [self.faidx.idx[x][0] for x in self.faidx.idx]
            if (self.type == 'DNA') or (self.type == 'RNA'):
                tmp_gc_at_gap = [(self.records[rname].GC(), self.records[rname].AT(), self.records[rname].gap()) 
                                 for rname in self.records]
                tmp_gc = [x[0] for x in tmp_gc_at_gap]
                tmp_at = [x[1] for x in tmp_gc_at_gap]
                tmp_gap = [x[2] for x in tmp_gc_at_gap]

            if each:
                # seq_len, %gc, %at, %gap for DNA/RNA
                # seq_len for Protein
                if (self.type == 'DNA') or (self.type == 'RNA'):
                    tmp_gc = ['%.2f' % (j/tmp_len[i]) for i, j in enumerate(tmp_gc)]
                    tmp_at = ['%.2f' % (j/tmp_len[i]) for i, j in enumerate(tmp_at)]
                    tmp_gap = ['%.2f' % (j/tmp_len[i]) for i, j in enumerate(tmp_gap)]
                    fout.write('seq\tseq_len\tGC\tAT\tGap\n')
                    for i, j in enumerate(self.faidx.idx.keys()):
                        fout.write('\t'.join(map(str, [j, tmp_len[i], tmp_gc[i], tmp_at[i], tmp_gap[i]])) + '\n')
                else:
                    fout.write('seq\tseq_len\n')
                    for i, j in enumerate(self.faidx.idx.keys()):
                        fout.write('\t'.join(map(str, [j, tmp_len[i]])) + '\n')
            else:
                # num_seqs, sum_len, %gc, min_len, max_len, avg_len, median_len, N50, N90, %gap
                num_seqs = self.len
                sum_len = sum(tmp_len)
                min_len = min(tmp_len)
                max_len = max(tmp_len)
                avg_len = mean(tmp_len)
                median_len = median(tmp_len)
                if (self.type == 'DNA') or (self.type == 'RNA'):
                    N50 = self.N50()
                    N90 = self.N50(0.9)
                    gc = '%.2f' % (sum(tmp_gc) / sum_len)
                    at = '%.2f' % (sum(tmp_at) / sum_len)
                    gap = '%.2f' % (sum(tmp_gap) / sum_len)
                    fout.write('file\ttype\tnum_seqs\tsum_len\tmin_len\tmax_len\tavg_len\tmedian_len\tN50\tN90\tGC\tAT\tgap\n')
                    fout.write('\t'.join(map(str, [self.filename, self.type, num_seqs, sum_len, min_len,
                                               max_len, avg_len, median_len, N50, N90, gc, at, gap])) + '\n')
                else:
                    fout.write('file\ttype\tnum_seqs\tsum_len\tmin_len\tmax_len\tavg_len\tmedian_len\n')
                    fout.write('\t'.join(map(str, [self.filename, self.type, num_seqs, sum_len, min_len,
                                                   max_len, avg_len, median_len])) + '\n')
            fout.close()


def prepare_argparser():
    """Prepare optparser object. New options will be added in this function first.
    """

    description = "%(prog)s -- a Python toolkit for fasta manipulation."
    epilog = "For command line options of each command, type: %(prog)s COMMAND -h"

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        description=description, epilog=epilog)
    subparsers = argparser.add_subparsers(dest='subcommand_name')

    # transform sequences (reverse, complement, extract ID...)
    add_seq_parser(subparsers)

    # get subsequences by region/gtf/bed, including flanking sequences
    add_subseq_parser(subparsers)

    # simple statistics of FASTA/Q files
    add_stats_parser(subparsers)

    return argparser

def add_out_option(parser):
    """add out options
    """
    group_output = parser.add_argument_group("Output files arguments")
    group_output.add_argument("-o", "--output", type=str, metavar="", default='-',
                        help='out file ("-" for stdout, suffix .gz for gzipped out)')

def add_input_option(parser):
    '''add input options
    '''
    group_input = parser.add_mutually_exclusive_group(required=True)
    group_input.add_argument('-fa', '--fasta', type=str, metavar='',
                        help='Fasta file input (gzipped file is supported)')
    group_input.add_argument('-fq', '--fastq', type=str, metavar='',
                        help='Fastq file input (gzipped file is supported)')

def add_fa_input_option(parser):
    '''add input for fasta
    '''
    group_input = parser.add_argument_group("Input files arguments")
    group_input.add_argument('-fa', '--fasta', type=str, metavar='', required=True,
                        help='Fasta file input (gzipped file is supported)')

def add_fq_input_option(parser):
    '''add input for fastq
    '''
    group_input = parser.add_argument_group("Input files arguments")
    group_input.add_argument('-fq', '--fastq', type=str, metavar='', required=True,
                        help='Fastq file input (gzipped file is supported)')

def add_report_flag(parser):
    '''add common flags
    '''
    group1 = parser.add_mutually_exclusive_group()
    group1.add_argument('-l', '--lower', action='store_true',
                        help='print sequences in lower case')
    group1.add_argument('-u', '--upper', action='store_true',
                        help='print sequences in upper case')
    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument('-s', '--seq', action='store_true',
                        help='only print sequences')
    group2.add_argument('-n', '--name', action='store_true',
                        help='only print names')
    group2.add_argument('-q', '--quality', action='store_true',
                        help='only print qualities')
    group3 = parser.add_mutually_exclusive_group()
    group3.add_argument('-rv', '--reverse', action='store_true',
                        help='reverse sequence')
    group3.add_argument('-rc', '--revcomp', action='store_true',
                        help='reverse complement sequence')
    parser.add_argument('-M', '--maxlen', type=int, metavar='', default=-1,
                        help='only print sequences shorter than the maximum length')
    parser.add_argument('-m', '--minlen', type=int, metavar='', default=0,
                        help='only print sequences longer than the minimum length')

def add_subseq_option(parser):
    '''add options for subseq
    '''
    group_subseq = parser.add_mutually_exclusive_group(required=True)
    group_subseq.add_argument('-bed', type=str, metavar='', 
                           help='tab-delimited BED file')
    group_subseq.add_argument('-gtf', type=str, metavar='',
                           help='GTF file')
    parser.add_argument('-ss', '--strand', action='store_true',
                        help='Force strandedness.')
    parser.add_argument('--up', type=int, metavar='', default=0,
                        help='up stream length')
    parser.add_argument('--down', type=int, metavar='', default=0,
                        help='down stream length')

def add_seq_parser(subparsers):
    """Add main function 'seq' argument parsers.
    """
    argparser_seq = subparsers.add_parser("seq",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                          help="transform sequences (reverse, complement, extract ID...)")
    # group for input files
    add_input_option(argparser_seq)
    # group flags for report
    group_report = argparser_seq.add_argument_group('Flags')
    add_report_flag(group_report)
    # group for output files
    add_out_option(argparser_seq)

def run_seq(args):
    '''
    function for 'seq' subcommand
    '''
    if args.fasta:
        fa = Fasta(args.fasta)
        fa.report(args.output, args.minlen, args.maxlen, args.reverse, args.revcomp,
                  args.lower, args.upper, args.name, args.seq)
    else:
        fq = Fastq(args.fastq)

def add_subseq_parser(subparsers):
    '''Add main function 'subseq' argument parsers.
    '''
    argparser_subseq = subparsers.add_parser("subseq",
                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                            help="get subsequences by gtf/bed of FASTA")
    # group for input files
    add_fa_input_option(argparser_subseq)
    # group for regions to sub
    add_subseq_option(argparser_subseq)
    # group flags for report
    group_report = argparser_subseq.add_argument_group('Report flags')
    add_report_flag(group_report)
    # group for output files
    add_out_option(argparser_subseq)

def run_subseq(args):
    '''function for 'subseq' command, only support fasta
    '''
    pass

def add_stats_parser(subparsers):
    """Add main function 'stats' argument parsers.
    """
    argparser_stats = subparsers.add_parser("stats",
                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                            help="simple statistics of FASTA/Q files")
    # group for input files
    add_input_option(argparser_stats)
    # group for stats
    argparser_stats.add_argument('--each', action='store_true',
                                  help='print stats for each sequence')
    # group for output files
    add_out_option(argparser_stats)

def run_stats(args):
    '''function for 'stats' function
    '''
    if args.fasta:
        fa = Fasta(args.fasta)
        fa.stats(args.output, args.each)
    else:
        pass


def main():

    argparser = prepare_argparser()

    if len(sys.argv) < 2:
        argparser.print_help()
        sys.exit(1)

    args = argparser.parse_args()
    subcommand = args.subcommand_name

    if subcommand == "seq":
        sys.exit(run_seq(args))

    elif subcommand == 'subseq':
        sys.exit(run_subseq(args))

    elif subcommand == 'stats':
        sys.exit(run_stats(args))


if __name__ == '__main__':

    try:
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(1)
