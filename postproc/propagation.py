#!/usr/bin/env python
"""
Eventually, this is to be a helper utility for determining what
information must be gathered to perform uncertainty quantification for
various derived quantities of interest.

TODO Still very much a slow work in progress
"""
from __future__ import division, print_function
from sympy.parsing.sympy_parser import parse_expr
import collections
import fileinput
import sympy
import sys
import tempfile

def parser(filenames):
    r'''
    Parse the provided filenames (or stdin if empty) into a
    collections.OrderedDict of symbol -> SymPy expression entries.
    See doctests for an example of the accepted syntax.

    >>> with tempfile.NamedTemporaryFile() as f:
    ...     f.write("a=1       # Comments       \n")
    ...     f.write(" b =  a+1 # Reuse earlier  \n")
    ...     f.write("c  = d+e  # Purely symbolic\n")
    ...     f.write("   f      # Nameless result\n")
    ...     f.flush()
    ...     parser(f.name)
    OrderedDict([('a', 1), ('b', 2), ('c', d + e), ('line4', f)])
    '''
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


def partials(e):
    r'''
    Given a SymPy expression e or string parsable as such, produce a
    collections.defaultdict where keys are e.free_symbols and values
    are the simplified derivative of e with respect to the key.
    The defaultdict returns sympy.Integer(0) for any key not in
    e.free_symbols.

    >>> partials("1")
    defaultdict(<class 'sympy.core.numbers.Integer'>, {})

    >>> partials("b**2 + a + 1")
    defaultdict(<class 'sympy.core.numbers.Integer'>, {b: 2*b, a: 1})
    '''
    if isinstance(e, basestring):
        e = parse_expr(e)
    r = collections.defaultdict(sympy.Integer)
    for s in e.free_symbols:
        r[s] = e.diff(s).simplify()
    return r

# def main(args):
#     symbol_table = parser([])
#     return 0
#
# if __name__=='__main__':
#     sys.exit(main(*sys.argv[1:]))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
