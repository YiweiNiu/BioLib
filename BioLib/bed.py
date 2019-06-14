#!/usr/bin/env python

'''
for python3
'''

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

import sys
import os
import argparse
from copy import deepcopy

# self module
from misc import cmp3

class Bed3:
    '''
    a class for bed3
    '''
    __slots__ = ('__n', '__Format', 'description', 'items', 'padding', '__chrom', '__start', '__end')

    def __init__(self, x, strict=True, n=3, description=None):
        '''
        init for all bed: n=3, bed3; n=6, bed6; n=12, bed12
        @parameter x: str of bed format, list, dict, tuple, or bed3 instance
        @parameter n: the type of bed
        @parameter strict: if strict, all fields should follow the types of bed
                           else, specific fields should follow
                           e.g. chrom, start, end, strand, block*
        '''
        self.__n = n
        if self.__n == 3:
            self.__Format = (str, int, int)
        elif self.__n == 6:
            self.__Format = (str, int, int, str, int, str)
        elif self.__n == 12:
            self.__Format = (str, int, int, str, int, str, int, int, \
                             int, int, tuple, tuple)
        else:
            raise Exception('Only 3/6/12 are supported value for bed instance.')

        if isinstance(x, str):
            lst = x.strip('\n').split()
            if len(lst) < 3:
                raise Exception('Three fields at least are needed to init bed class.')
        elif isinstance(x, list):
            if len(x) < 3:
                raise Exception('Three fields at least are needed to init bed class.')
            lst = deepcopy(x)
        elif isinstance(x, dict):
            if len(x) < 3:
                raise Exception('Three fields at least are needed to init bed class.')
            lst = [x[i] for i in x]
        elif isinstance(x, tuple):
            if len(x) < 3:
                raise ValueError('Three fields at least are needed to init bed class.')
            lst = list(x)
        elif isinstance(x, type(self)):
            lst = [getattr(x, i) for i in self.__slots__]
        else:
            raise ValueError('%s initilization needs string, list, dict, tuple or bed3 class.')

        items = []
        for i in range(self.__n):
            if strict:
                try:
                    if i in (10, 11):
                        items.append(self.__Format[i](map(int, (i for i in lst[i].split(',') if i!=''))))
                    else:
                        items.append(self.__Format[i](lst[i]))
                except:
                    raise Exception('Bed%s instance initilization error.' %self.__n)
            else:
                try:
                    if i == 0:
                        items.append(str(lst[i]))
                    elif i == 1:
                        items.append(int(lst[i]))
                    elif i == 2:
                        items.append(int(lst[i]))
                    elif i == 5:    # strand
                        items.append(str(lst[i]))
                    elif i == 9:    # blockCount
                        items.append(int(lst[i]))
                    elif i in (10, 11):
                        items.append(tuple(map(int, (i for i in lst[i].split(',') if i!=''))))
                    else:
                        items.append(lst[i])
                except:
                    raise Exception('Bed%s instance initilization error.' %self.__n)

        if items[2] <= items[1]:
            raise Exception('Interval end %s is larger than interval start %s.' %(items[2], items[1]))

        # other fields that do not include in the class
        items.append('\t'.join(map(str, lst[self.__n:])))

        self.items = items
        self.chrom = self.items[0]
        self.start = self.items[1]
        self.end = self.items[2]
        self.padding = self.items[-1]   # other fields
        self.description = description if description else ''

    @property
    def chrom(self):
        return self.__chrom
    @chrom.setter
    def chrom(self, value):
        assert isinstance(value, str)
        self.__chrom = value
        self.items[0] = value

    @property
    def start(self):
        return self.__start
    @start.setter
    def start(self, value):
        assert isinstance(value, int)
        self.__start = value
        self.items[1] = value

    @property
    def end(self):
        return self.__end
    @end.setter
    def end(self, value):
        assert isinstance(value, int)
        self.__end = value
        self.items[2] = value

    def __str__(self):  # bed3 and bed6
        return '\t'.join(map(str, self.items))

    def __repr__(self):
        return 'Bed'+str(self.__n)+' object:\n'+str(self)+'\n'

    def __len__(self): #All bed
        return self.end - self.start

    def __eq__(self, x):
        '''self == x
        '''
        if not isinstance(x, type(self)):
            raise TypeError('Object provided is not a Bed instance.')
        return (self.chrom == x.chrom) and (self.start == x.start) and (self.end == x.end)

    def __ne__(self, x):
        '''self != x
        '''
        if not isinstance(x, (Bed3, Bed6, Bed12)):
            raise TypeError('Object provided is not a Bed instance.')
        return (self.chrom != x.chrom) or (self.start != x.start) or (self.end != x.end)

    def __lt__(self, x):
        '''self < x
        '''
        if not isinstance(x, (Bed3, Bed6, Bed12)):
            raise TypeError('Object provided is not a Bed instance.')
        return (self.chrom < x.chrom) or (self.start < x.start) or (self.end < x.end)

    def __le__(self, x):
        '''self <= x
        '''
        if not isinstance(x, (Bed3, Bed6, Bed12)):
            raise TypeError('Object provided is not a Bed instance.')
        return self.__lt__(x) or self.__eq__(x)

    def __gt__(self, x):
        '''self > x
        '''
        if not isinstance(x, (Bed3, Bed6, Bed12)):
            raise TypeError('Object provided is not a Bed instance.')
        return (self.chrom > x.chrom) or (self.start > x.start) or (self.end > x.end)

    def __ge__(self, x):
        '''self >= x
        '''
        if not isinstance(x, (Bed3, Bed6, Bed12)):
            raise TypeError('Object provided is not a Bed instance')
        return self.__gt__(x) or self.__eq__(x)

    def expand(self, upstream=0, downstream=0):
        '''expand bed regions
        '''
        new_bed = deepcopy(self)
        if hasattr(new_bed, 'strand'):
            if (new_bed.strand == '.') or (new_bed.strand == '+'):
                new_bed.start -= upstream if (new_bed.start-upstream > 0) else 0
                new_bed.end += downstream
            else:
                new_bed.end += upstream
                new_bed.start -= downstream if (new_bed.start-downstream > 0) else 0
        else:
            new_bed.start -= upstream if (new_bed.start-upstream) > 0 else 0
            new_bed.end += downstream

    def __is_overlap(self, x, split=False):
        '''helper function for is_overlap
        '''
        if split:
            if (not hasattr(self, 'blockSizes')) or (not hasattr(self, 'blockStarts')) \
               or (not hasattr(x, 'blockSizes')) or (not hasattr(x, 'blockStarts')):
                logger.warning('Object does not contain block information, and split parameter will be ignored.')
                return self.__is_overlap(x, False)
            elif (not hasattr(self, 'blockCount')) or (not hasattr(self, 'blockCount')):
                logger.warning('Object does not contain block information, and split parameter will be ignored.')
                return self.__is_overlap(x, False)
            elif self.blockCount==1 and x.blockCount==1:
                return self.__is_overlap(x, False)
            else:
                tmp_a, tmp_b = list(), list()
                for i, j in enumerate(self.blockSizes):
                    tmp_a.append(Bed3((self.chrom, self.start+self.blockStarts[i], self.start+self.blockStarts[i]+j)))
                for i, j in enumerate(x.blockSizes):
                    tmp_b.append(Bed3((x.chrom, x.start+x.blockStarts[i], x.start+x.blockStarts[i]+j)))
                for i in tmp_a:
                    for j in tmp_b:
                        if i.__is_overlap(j, False):
                            return True
                return False
        else:
            if set(range(self.start, self.end)).intersection(set(range(x.start, x.end))):
                return True
            else:
                return False

    # for all bed
    def is_overlap(self, x=None, strandedness=False, split=False):
        '''
        see if two region overlap with each other
        '''
        if not isinstance(x, (Bed3, Bed6, Bed12)):
            raise TypeError('Object provided is not a %s instance.' %self.__class__)

        if self.chrom != x.chrom:
            return False

        if strandedness:
            if not hasattr(self, 'strand'):
                raise Exception('Object does not have strand attribute.')
            elif not hasattr(x, 'strand'):
                raise Exception('Object does not have strand attribute.')
            elif self.strand=='.' or x.strand=='.':
                logger.warning('Strand information will be ignored, since strand is empty.')
                return self.is_overlap(x, False, split)
            elif self.strand != x.strand:
                return False
            else:
                if split:
                    return self.__is_overlap(x, True)
                else:
                    return self.__is_overlap(x, False)
        else:
            if split:
                return self.__is_overlap(x, True)
            else:
                return self.__is_overlap(x, False)

    def distance(self, x, strandedness=False):
        '''
        get the distance of self and x
        used in merge operation
        '''
        if self.chrom != x.chrom:
            return None

        if self.is_overlap(x, strandedness):
            return 0

        if strandedness:
            if not hasattr(self, 'strand'):
                raise Exception('Object does not have strand attribute.')
            elif not hasattr(x, 'strand'):
                raise Exception('Object does not have strand attribute.')
            elif self.strand=='.' or x.strand=='.':
                logger.warning('Strand information will be ignored, since strand is empty.')
                return self.distance(x, False)
            elif self.strand != x.strand:
                return None
            else:
                tmp_list = sorted([self.start, self.end, x.start, x.end])
                return tmp_list[2]-tmp_list[1]
        else:
            tmp_list = sorted([self.start, self.end, x.start, x.end])
            return tmp_list[2]-tmp_list[1]


class Bed6(Bed3):
    '''
    a class for bed6
    '''
    __slots__ = ('__name', '__score', '__strand')

    def __init__(self, x, strict=True, n=6, description=None):
        Bed3.__init__(self, x, strict, n, description=None)
        self.name = self.items[3]
        self.score = self.items[4]
        self.strand = self.items[5]

    @property
    def name(self):
        return self.__name
    @name.setter
    def name(self, value):
        self.__name = value
        self.items[3] = value

    @property
    def score(self):
        return self.__score
    @score.setter
    def score(self, value):
        self.__score = value
        self.items[4] = value

    @property
    def strand(self):
        return self.__strand
    @strand.setter
    def strand(self, value):
        if value not in ('+', '-', '.'):
            raise ValueError('Strand must be ., +, or -.')
        self.__strand = value
        self.items[5] = value


class Bed12(Bed6):
    '''
    a class for bed12
    '''
    __slots__ = ('__thickStart', '__thickEnd', '__itemRgb', '__blockCount', '__blockSizes', '__blockStarts')

    def __init__(self, x, strict=True, n=12, description=None):
        Bed6.__init__(self, x, strict, n=12, description=None)
        self.thickStart = self.items[6]
        self.thickEnd = self.items[7]
        self.itemRgb = self.items[8]

        if self.items[9] != len(self.items[10]) != len(self.items[11]):
            raise Exception('Bed%s instance initilization error. Block information error.' %n)

        self.blockCount = self.items[9]
        self.blockSizes = self.items[10]
        self.blockStarts = self.items[11]

    @property
    def thickStart(self):
        return self.__thickStart
    @thickStart.setter
    def thickStart(self, value):
        self.__thickStart = value
        self.items[6] = value

    @property
    def thickEnd(self):
        return self.__thickEnd
    @thickEnd.setter
    def thickEnd(self, value):
        self.__thickEnd = value
        self.items[7] = value

    @property
    def itemRgb(self):
        return self.__itemRgb
    @itemRgb.setter
    def itemRgb(self, value):
        self.__itemRgb = value
        self.items[8] = value

    @property
    def blockCount(self):
        return self.__blockCount
    @blockCount.setter
    def blockCount(self, value):
        assert isinstance(value, int)
        self.__blockCount = value
        self.items[9] = value

    @property
    def blockSizes(self):
        return self.__blockSizes
    @blockSizes.setter
    def blockSizes(self, value):
        self.__blockSizes = value
        self.items[10] = value

    @property
    def blockStarts(self):
        return self.__blockStarts
    @blockStarts.setter
    def blockStarts(self, value):
        self.__blockStarts = value
        self.items[11] = value

    def __str__(self):  # bed12
        return '\t'.join(list(map(str, self.items[:10])) + \
                [','.join(map(str, self.items[10])) + ','] + \
                [','.join(map(str, self.items[11])) + ','] + \
                [self.padding])


class BedList:
    '''class for list of bed instance
    '''
    __slots__= ['records']
    def __init__(self, bed_lst):
        '''
        @parameter bed_lst: a list or generator of bed instance
        '''
        self.records = []  # not check the validity here for now
        for record in bed_lst:
            self.records.append(record)

    def __iter__(self):
        for bed in self.records:
            yield bed

    def __len__(self):
        return len(self.records)

    def __add__(self, x):
        '''concat two BedList
        '''
        if not isinstance(x, type(self)):
            raise Exception('Input object %s is not a BedList instance.' %x)
        return self.__class__(self.records + x.records)

    def sort(self, reverse=False):
        '''sort the list
        '''
        self.records.sort(reverse=reverse)

    def merge(self, strandedness=False, distance=0, sorted=False):
        '''merge overlapped intervals
        '''
        if not sorted:
            self.sort()

        stack = list()
        stack.append(self.records[0])

        for i in range(1, len(self)):
            top = stack[-1]
            dist = top.distance(self.records[i], strandedness=strandedness)
            if dist == None:
                stack.append(self.records[i])
            elif dist <= distance:
                if top.end < self.records[i].end:
                    top.end = self.records[i].end
                    stack.pop()
                    stack.append(top)
            else:
                stack.append(self.records[i])

        return self.__class__(stack)

    def report(self, outfilename='-'):
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
            for bed in self.records:
                fout.write(str(bed) + '\n')
            fout.close()


class BedFile:
    '''a class for bed file, read bed record from file
    '''
    def __init__(self, filename, strict=False, bed_type=None, sorted=False):
        '''
        init a bed class

        @parameter finename: file name of the fasta file
        @parameter strict: whether to use strict mode to open file
        @parameter bed_type: the type of bed file (bed3, bed6, bed12)
        '''
        self.filename = os.path.abspath(os.path.realpath(filename))
        self.mode = strict  # mode to create bed instance

        if bed_type:
            if bed_type not in ('bed3', 'bed6', 'bed12'):
                raise Exception('Only bed3/bed6/bed12 are valid bed_type')
            self.bed_type = bed_type

        if filename.lower().endswith('.gz'):
            import gzip
            self._file_opener = gzip.open
        else:
            self._file_opener = open

        # check the bed file, get bed_type if it's None
        self.__check()

        if self.bed_type == 'bed3':
            self.__bed_class = Bed3
        elif self.bed_type == 'bed6':
            self.__bed_class = Bed6
        else:
            self.__bed_class = Bed12

        self.records = BedList(self.__to_record())

    def __check(self):
        '''use the minimum number of fields to judge the type of bed to use
        '''
        try:
            fin = self._file_opener(self.filename, 'rt')
        except IOError:
            raise Exception('Error in opening %s' % self.filename)
        else:
            __field_num = set()
            for line in fin:
                line = line.strip()
                if line.startswith('#'):
                    continue
                line = line.split()
                if len(line) < 3:
                    raise Exception('\t'.join(line) + ' does not have three fields.')
                try:
                    start, end = int(line[1]), int(line[2])
                except:
                    raise Exception('The start %s and end %s should be int.' %(line[1], line[2]))

                if start >= end:
                    raise Exception('The start %s should be less than the end %s.' %(line[1], line[2]))

                __field_num.add(len(line))

            if len(__field_num) >= 2:
                logger.warning('The file %s has variable lengths.' %self.filename)

            if hasattr(self, 'bed_type'):
                if self.bed_type == 'bed6':
                    if min(__field_num) < 6:
                        raise Exception('Some lines of %s do not have 6 fields.' %self.filename)
                elif self.bed_type == 'bed12':
                    if min(__field_num) < 12:
                        raise Exception('Some lines of %s do not have 12 fields.' %self.filename)
            else:
                if min(__field_num) >= 12:
                    self.bed_type = 'bed12'
                elif min(__field_num) >= 6:
                    self.bed_type = 'bed6'
                else:
                    self.bed_type = 'bed3'

    def __to_record(self):
        '''
        convert fasta file to record
        modified from: https://www.biostars.org/p/710/#120760
        '''
        try:
            fin = self._file_opener(self.filename, 'rt')
        except IOError:
            raise Exception('Error in opening %s' % self.filename)
        else:
            for line in fin:
                line = line.strip()
                if line.startswith('#'):
                    continue
                yield self.__bed_class(line, self.mode)

    def __len__(self):
        return len(self.records)

    def __repr__(self):
        return '%s("%s")' % (self.__bed_class, self.filename)

    def __iter__(self):
        for bed in self.records:
            yield bed

    def sort(self, reverse=False):
        '''
        '''
        self.records.sort(reverse=reverse)

    def merge(self, x=None, strandedness=False, distance=0, sorted=False):
        '''merge
        '''
        if x:
            if not isinstance(x, type(self)):
                raise Exception('The object to merge is not a BedFile instance.')
            __merged_interval = self.records + x.records
            merged_interval = __merged_interval.merge(strandedness, distance, sorted)
        else:
            merged_interval = self.records.merge(strandedness, distance, sorted)
        return BedList(merged_interval)


def prepare_argparser():
    """Prepare optparser object. New options will be added in this function first.
    """

    description = "%(prog)s -- a Python toolkit for bed manipulation."
    epilog = "For command line options of each command, type: %(prog)s COMMAND -h"

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        description=description, epilog=epilog)
    subparsers = argparser.add_subparsers(dest='subcommand_name')

    # merge
    add_merge_parser(subparsers)

    return argparser

def add_out_option(parser):
    """add out options
    """
    group_output = parser.add_argument_group("Output files arguments")
    group_output.add_argument("-o", "--output", type=str, metavar="", default='-',
                        help='out file ("-" for stdout, suffix .gz for gzipped out)')

def add_single_input(parser):
    '''add single bed input options
    '''
    group_input = parser.add_argument_group("Input files arguments")
    group_input.add_argument('-i', type=str, metavar='', required=True, dest='input',
                             help='Bed file input (gzipped file is supported)')

def add_merge_parser(subparsers):
    """Add main function 'merge' argument parser.
    """
    argparser_merge = subparsers.add_parser("merge",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                          help="transform sequences (reverse, complement, extract ID...)")
    # group for input files
    add_single_input(argparser_merge)
    # group flags for mrege
    group_mrege = argparser_merge.add_argument_group('Merge options')
    group_mrege.add_argument('-d', type=int, metavar='', default=0, dest='distance',
        help='maximum distance between features allowed for features to be merged.')
    group_mrege.add_argument('-s', action='store_true', dest='strandedness',
        help='force strandedness.')
    # group for output files
    add_out_option(argparser_merge)


def run_merge(args):
    '''
    function for 'merge' subcommand
    '''
    bed_file = BedFile(args.input)
    merged_interval = bed_file.merge(None, args.strandedness, args.distance, False)
    merged_interval.report(args.output)


def test():
    '''
    test function
    '''
    # bed3
    a3 = Bed3('chr1\t1\t10\ta3')
    b3 = Bed3('chr1\t6\t15\tb3')
    # bed6
    a1 = Bed6('chr1\t100\t200\ta1\t1\t+')
    a2 = Bed6('chr1\t180\t250\ta2\t2\t+')
    print(a1.chrom, a2.chrom)
    a1.chrom = 'chrM'
    print(a1.chrom)
    print(a1)

def main():

    argparser = prepare_argparser()

    if len(sys.argv) < 2:
        argparser.print_help()
        sys.exit(1)

    args = argparser.parse_args()
    subcommand = args.subcommand_name

    if subcommand == "merge":
        sys.exit(run_merge(args))

    elif subcommand == 'subseq':
        sys.exit(run_subseq(args))

    elif subcommand == 'stats':
        sys.exit(run_stats(args))


if __name__ == '__main__':

    try:
        #test()
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(1)


