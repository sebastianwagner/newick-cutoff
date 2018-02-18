#!/usr/bin/env python3

# -*- coding: utf-8 -*-

from Bio import Phylo
from Bio.Phylo import NewickIO  # directly called due to stdout usage
import sys
import logging

def perclade(clade):
    return clade

def walktree(clade):
    if clade.name:
        clade = perclade(clade)
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
    if treefile == '-':
        return NewickIO.parse(sys.stdin)
    return Phylo.parse(treefile, 'newick')


def main():
    stdout_handler = logging.StreamHandler(sys.stderr)
    handlers = [stdout_handler]
    logging.basicConfig(
        level=logging.INFO,
        format='[%(levelname)s - %(message)s]',
        handlers=handlers
    )
    log = logging.getLogger('LOGGER_NAME')

    treefile = sys.argv[1]
    trees = readtrees(treefile)

    trees = relabeltree(trees)

    NewickIO.write(trees, sys.stdout)


if __name__ == '__main__':
    main()

