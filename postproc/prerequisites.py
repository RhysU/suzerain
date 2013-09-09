#!/usr/bin/env python
"""
Eventually, this is to be a helper utility for determining
what information must be gathered to perform uncertainty
quantification for various derived quantities of interest.
"""
from __future__ import division, print_function
from sympy.parsing.sympy_parser import parse_expr
import fileinput
import collections
import sys

def parser(filenames):
    """
    Parse the provided filenames (or stdin if empty) into
    an OrderedDict of symbol -> SymPy expression entries.
    """
    # Accumulate symbol definitions maintaining declaration order
    symbol_table = collections.OrderedDict()

    # Process input line-by-line...
    try:
        for line in fileinput.input(filenames):

            # ...remove comments which occur after the first '#' character
            head, sep, tail = line.partition('#')
            line = head if head else tail

            # ...split lines into "token = rhs" or just "rhs" with the
            # latter implicitly naming the result like "line1234"
            symbol, sep, expr = line.partition('=')
            expr = expr.strip()
            if not expr:
                symbol, expr = 'line' + `fileinput.lineno()`, symbol
            symbol = symbol.strip()
            if symbol in symbol_table:
                raise LookupError("Duplicate symbol '%s' detected at %s:%d" % (
                                symbol, fileinput.filename(),
                                fileinput.filelineno()))

            # ...parse and save the result into the known symbol table
            # augmenting any parsing errors with location information
            if expr:
                try:
                    symbol_table[symbol] = parse_expr(expr, symbol_table)
                except SyntaxError as e:
                    e.filename = fileinput.filename()
                    e.lineno   = fileinput.lineno()
                    raise

    # ...being sure to clean up after our use of fileinput
    finally:
        fileinput.close()

    return symbol_table

# def main(args):
#     return 0
#
# if __name__=='__main__':
#     sys.exit(main(*sys.argv[1:]))
