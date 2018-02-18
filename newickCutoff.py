#!/usr/bin/env python3

# -*- coding: utf-8 -*-

from Bio import Phylo
from Bio.Phylo import NewickIO  # directly called due to stdout usage
import sys
import logging
import argparse

cutoff = 75

def perclade(clade):
    if clade.confidence:
        if clade.confidence < cutoff:
            clade.confidence = None
    return clade

def walktree(clade):
    clade = perclade(clade)
    if clade.name:
        pass
    # fix broken confidence writer when working with PAUP* style files.
    # TreeGraph could not handle confidence behind colons
    # @link http://wiki.christophchamp.com/index.php?title=Newick_phylogenetic_tree_format
    # @link https://en.wikipedia.org/wiki/Newick_format
    # @see Bio.Phylo.NewickIO.Writer#_info_factory
    elif clade.confidence is not None:
        clade.name = str(clade.confidence)
        clade.confidence = None
    if clade.clades:
        for subclade in clade.clades:
            walktree(subclade)

def relabeltree(trees):
    outtrees = []
    for tree in trees:
        if tree.clade:
            walktree(tree.clade)
        outtrees.append(tree)
    return outtrees

def readtrees(treefile):
    return Phylo.parse(treefile, 'newick')


def main():
    prog = sys.argv[0]
    description = ('Parse newick tree an perform action on'
                   'each non-root node')
    parser = argparse.ArgumentParser(prog=prog, description=description)
    parser.add_argument('infile', nargs='?', type=argparse.FileType(),
                        help='a Newick treefile')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        help='changed Newick outfile')
    parser.add_argument('--cutoff', dest='cutoff', nargs='?', type=int, default=75,
                        help='value at or beneath which no inner node'
                             'confidences are snown any more')
    options = parser.parse_args()

    stdout_handler = logging.StreamHandler(sys.stderr)
    handlers = [stdout_handler]
    logging.basicConfig(
        level=logging.INFO,
        format='[%(levelname)s - %(message)s]',
        handlers=handlers
    )
    log = logging.getLogger('LOGGER_NAME')

    infile = options.infile or sys.stdin
    outfile = options.outfile or sys.stdout
    cutoff = options.cutoff
    trees = readtrees(infile)
    trees = relabeltree(trees)

    NewickIO.write(trees, outfile)


if __name__ == '__main__':
    main()

