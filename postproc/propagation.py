#!/usr/bin/env python
r'''
Eventually, this is to be a helper module for determining what information
must be gathered to perform uncertainty quantification for various
derived quantities of interest.

The underlying model is
  \vec{d} = \vec{x} - \vec{\beta} - \vec{\epsilon}
where
  \vec{x} is some deterministic truth
  \vec{d} is one observation of \vec{x}
  \vec{\epsilon} is a zero-mean, normally-distributed measurement error
  \vec{\beta} is a bias error which is small relative to \vec{\epsilon}
Additionally, \vec{\beta} is assumed to be independent of \vec{\epsilon}.
Assume also \vec{\epsilon} is a zero-mean, normally-distributed random
variable with some known covariance matrix \Sigma containing scalar
components \sigma_{ij}.

TODO Still very much a slow work in progress
'''
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


def partials(f):
    r'''
    Given a SymPy expression f or any string parsable as such, produce
    a defaultdict df where referencing df[x] produces the precomputed
    result f.diff(x).simplify() for any x.

    >>> a, b, c = sympy.symbols('a, b, c')
    >>> df = partials(b**2 + a + 1)
    >>> df[a], df[b], df[c]
    (1, 2*b, 0)

    >>> df = partials("1 + 2 + 3")
    >>> df.keys()
    []
    '''
    if isinstance(f, basestring):
        f = parse_expr(f)
    df = collections.defaultdict(lambda: sympy.Integer(0))
    for x in f.free_symbols:
        df[x] = f.diff(x).simplify()
    return df


def mixed_partials(f, df=None):
    r'''
    Given a SymPy expression f or any string parsable as such, produce a
    defaultdict of defaultdicts ddf where referencing ddf[x][y] produces
    the precomputed result f.diff(x,y).simplify() for any x and y.

    >>> a, b, c = sympy.symbols('a, b, c')
    >>> ddf = mixed_partials(a**2 + a*b + b**2 + 1)
    >>> ddf[a][a], ddf[a][b], ddf[b][a], ddf[b][b]
    (2, 1, 1, 2)
    >>> ddf[a][c], ddf[c][a], ddf[c][c]
    (0, 0, 0)
    '''
    ddf = collections.defaultdict(
            lambda: collections.defaultdict(lambda: sympy.Integer(0))
        )
    if df is None:
        df = partials(f)
    for (x, dfdx) in df.iteritems():
        ddf[x] = partials(dfdx)

    return ddf

def prerequisites(f, df=None, ddf=None):
    r'''
    Given a SymPy expression f or any string parsable as such, return
    a tuple (E, Cov) where E and Cov are is the set of free symbol
    expectations and covariances, respectively, necessary to compute an
    estimate of E[f] and Var[f] using Taylor Series Methods (TSM).

    Applying TSM to the underlying model yields

        E[f(x)] = f(d) + (1/2) \sum_{i,j} \sigma_{ij} f_{i,j}(d)

    '''
    if isinstance(f, basestring):
        f = parse_expr(f)
    if df is None:
        df = partials(f)
    if ddf is None:
        ddf = mixed_partials(f, df)
    E, Cov = set(), set()

    return E, Cov

# def main(args):
#     symbol_table = parser([])
#     return 0
#
# if __name__=='__main__':
#     sys.exit(main(*sys.argv[1:]))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
