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
import itertools
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

# See section "Estimating uncertainty in derived quantities" in Suzerain's
# writeups/perfectgas.tex for derivation of the TSM results appearing below.
def prerequisites(f, df=None, ddf=None):
    r'''
    Given a SymPy expression f or any string parsable as such, return a
    list wherein unique tuples represents moments necessary for computing
    an estimate of E[f] and Var[f] using Taylor Series Methods (TSM).

    Applying TSM to the underlying model yields

        E[f(x)]   ~= f(d) + (1/2) \sum_{i,j} \sigma_{ij} f_{,ij}(d)

        Var[f(x)] ~=        \sum_{ i } \sigma_{ii} f_{,i}^2(d)
                   +    2   \sum_{i<j} \sigma_{ij} f_{,i}(d) f_{,j}(d)
                   -  f(d)  \sum_{i,j} \sigma_{ij} f_{,ij}(d)
                   - (1/4) (\sum_{i,j} \sigma_{ij} f_{,ij}(d))^2

    where f_{,i} denotes partial differentiation with respect to scalar
    component x_i and f_{,ij} denotes differentiation with respect to
    both components x_i and x_j.

    >>> prerequisites('1')
    []

    >>> prerequisites('log(x)')
    [(x,), (x, x)]

    >>> prerequisites('2*x + 3*y')
    [(x,), (x, x), (x, y), (y,), (y, y)]

    >>> prerequisites('x*y')
    [(x,), (x, x), (x, y), (y,), (y, y)]
    '''
    if isinstance(f, basestring):
        f = parse_expr(f)
    f = f.simplify()
    if df is None:
        df = partials(f)
    if ddf is None:
        ddf = mixed_partials(f, df)

    # Implementation heavily relies on set addition semantics combined
    # with the fact that all derivatives have been precomputed prior
    # to iteration.  Because the caller might know something we do not
    # about he or she wants to compute a subexpression (hinted via
    # df or ddf arguments), we look at all possible terms rather than
    # removing those which can be eliminated by smoothness or symmetry.
    m = set()

    # Quantities necessary to compute E[f(x)]
    ## Term:    f(d)
    for s in f.free_symbols:
        m.add((s,))
    ## Term: +  (1/2) \sum_{i,j} \sigma_{ij} f_{,ij}(d)
    for i in ddf.keys():
        for j in ddf[i].keys():
            f_ij = ddf[i][j]
            if not f_ij.is_zero:
                for s in f_ij.free_symbols:
                    m.add((s,))
                m.add(tuple(sorted([i, j]))) # Canonicalize

    # Quantities additionally necessary to compute Var[f(x)]
    ## Term:          \sum_{ i } \sigma_{ii} f_{,i}^2(d)
    for i in df.keys():
        f_i = df[i]
        if not f_i.is_zero:
            for s in f_i.free_symbols:
                m.add((s,))
            m.add((i, i))                    # Canonical
    ## Term: +    2   \sum_{i<j} \sigma_{ij} f_{,i}(d) f_{,j}(d)
    for (i, j) in itertools.combinations(df.keys(), 2):
        fifj = (df[i] * df[j]).simplify()
        if not fifj.is_zero:
            for s in fifj.free_symbols:
                m.add((s,))
            m.add(tuple(sorted([i, j])))     # Canonicalize
    ## These terms require no additional data relative to E[f(x)]:
    ## Term: -  f(d)  \sum_{i,j} \sigma_{ij} f_{,ij}(d)
    ## Term: - (1/4) (\sum_{i,j} \sigma_{ij} f_{,ij}(d))^2

    return sorted(m)

# TODO Implement something returning simplified expressions
# Test cases should be based on Table 1 of the article H. H. Ku. Notes
# on the use of propagation of error formulas. Journal of Research
# of the National Bureau of Standards. Section C: Engineering and
# Instrumentation, 70C(4):263-273, October 1966.  ISSN 0022-4316.

# def main(args):
#     symbol_table = parser([])
#     return 0
#
# if __name__=='__main__':
#     sys.exit(main(*sys.argv[1:]))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
