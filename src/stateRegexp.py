#!/usr/bin/env python3

import argparse
import contextlib
import sys
from array import array
from typing import TextIO

import regex as re


@contextlib.contextmanager
def smart_open(filename: str = None) -> TextIO:
    # along lines of https://stackoverflow.com/questions/17602878/
    if filename and filename != '-':
        fh = open(filename, 'wt')
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


class StateRegexp(object):
    BIVALENT = r"[KLNOPQR]([KLNOPQR]G+[KLNOPQR])[KLNOPQR]"
    REPRPC = r"[MNOPQR][KLNOPQR]([KLNOPQR][LK]+[KLNOPQR])[KLNOPQR][MNOPQR]"

    # input files are expected to have file names:
    #   ...<tissue>_<timepoint>_mm10_...
    # e.g. /some/path/lung_15.5_mm10_18_posterior.bed

    # input files should be basically BED with posterior probability
    # from chromhmm in the fifth column, e.g.:
    #   chr10   1000    1200    E4  0.991

    def __init__(self):
        self.alphabet = {
            "E16": "A",  # Tss active_tss E16
            "E18": "B",  # TssFlnk flanking_active_tss E18
            "E10": "C",  # Tx transcription E10
            "E8":  "D",  # TxWk weak_transcription E8
            "E14": "E",  # Enh enhancer E14
            "E13": "F",  # EnhLo weak_enhancer E13
            "E17": "G",  # TssBiv bivalent_tss E17
            "E12": "H",  # EnhPois poised_enhancer E12
            "E15": "I",  # EnhPr primed_enhancer E15
            "E11": "J",  # EnhG genic_enhancer E11
            "E1":  "K",  # ReprPC polycomb_repressed E1
            "E2":  "L",  # ReprPCWk polycomb_repressed_weak E2
            "E9":  "M",  # Het heterochromatin E9
            "E5":  "N",  # QuiesG quiescent_gene E5
            "E3":  "O",  # Quies quiescent E3
            "E4":  "P",  # Quies2 quiescent2 E4
            "E6":  "Q",  # Quies3 quiescent3 E6
            "E7":  "R",  # Quies4 quiescent4 E7
        }
        self.threshold = 0.50

    def _readfile(self, fnp, size=13627678) -> array:
        # read file into a character array (on char for every window)
        s = array('b', '.'.encode('ascii') * size)  # initialize w/ dots
        with open(fnp, "r") as f:
            i = 0
            for line in f:
                line = line.rstrip().split('\t')
                if float(line[4]) > 0.5:  # if p > 0.5, use alphabet
                    s[i] = ord(self.alphabet[line[3]])
                i += 1
        return s

    @staticmethod
    def mouse_tissue(fnp: str) -> str:
        # get tissue from filename
        return re.sub(r"^.*/(.*)_([0-9.]+|unknown)_mm10.*$", "\\1", fnp)

    @staticmethod
    def mouse_timepoint(fnp: str) -> str:
        # get time point from filename
        return re.sub(r"^.*/(.*)_([0-9.]+|unknown)_mm10.*$", "\\2", fnp)

    def regex(self, fnp: str, pattern: str, slop: int = 0, binsize: int = 200, tag: str = None, overlapped: bool = True,
              outfnp: str = None, bed3: bool = False):
        tissue = self.mouse_tissue(fnp)
        timepoint = self.mouse_timepoint(fnp)
        with smart_open(outfnp) as of:
            with open(fnp, 'rt') as f:  # open file
                n = 0
                s = self._readfile(fnp).tostring().decode('ascii')
                for i in re.finditer(pattern, s, overlapped=overlapped):  # find all positions matching regex
                    for reg in i.regs[1:]:
                        length = reg[1] - reg[0]
                        start = reg[0]
                        while n <= start:
                            line = f.readline()  # seek to right position
                            n += 1
                        line = line.rstrip().split('\t')[0:3]
                        line[2] = int(line[2]) + binsize * (length - 1)
                        if slop:
                            line[1] = max(0, int(line[1]) - slop)
                            line[2] = int(line[2]) + slop
                        if bed3:
                            of.write('\t'.join(map(str, line)) + '\n')
                        else:
                            if tag:
                                of.write(
                                    '\t'.join(map(str, line)) + '\t' + tissue + '\t' + timepoint + '\t' + tag + '\n')
                            else:
                                of.write('\t'.join(map(str, line)) + '\t' + tissue + '\t' + timepoint + '\n')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--type", "-t", choices=["bivalent", "reprpc"], required=True)
    parser.add_argument("--file", "-f", type=str, required=True)
    parser.add_argument("--out", "-o", type=str)
    parser.add_argument("--bed3", "-3", action="store_true", default=False)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if "bivalent" == args.type:
        StateRegexp().regex(args.file, StateRegexp.BIVALENT, 200, tag="bivalent",
                            overlapped=True, outfnp=args.out, bed3=args.bed3)
    elif "reprpc" == args.type:
        StateRegexp().regex(args.file, StateRegexp.REPRPC, 200, tag="reprpc",
                            overlapped=True, outfnp=args.out, bed3=args.bed3)
