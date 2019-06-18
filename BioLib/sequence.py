#!/usr/bin/env python

'''
for python3
'''

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

import sys
import os

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


