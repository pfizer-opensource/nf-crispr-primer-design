"""Utility functions."""

_DNA_complement = str.maketrans('ACGTacgt', 'TGCAtgca')

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_DNA_complement)[::-1]
