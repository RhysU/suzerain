#!/usr/bin/env python
r'''
Symbolic routines for obtaining uncertainty propagation expressions
based upon Taylor Series Methods (TSM).

The underlying model is

  \vec{d} = \vec{x} - \vec{\beta} - \vec{\epsilon}

where

  \vec{x} is some deterministic truth
  \vec{d} is one observation of \vec{x}
  \vec{\epsilon} is a zero-mean, normally-distributed measurement error
  \vec{\beta} is a bias error which is small relative to \vec{\epsilon}

Additionally, \vec{\beta} is assumed to be independent of \vec{\epsilon}.
Assume also \vec{\epsilon} is a zero-mean, normally-distributed random variable
with some known covariance matrix \Sigma containing scalar components
\sigma_{ij}.

Applying TSM to the underlying model yields

    E[f(x)]   ~=  f(d)
               + (1/2) \sum_{ i } \sigma_{ii} f_{,ii}(d)
               +       \sum_{i<j} \sigma_{ij} f_{,ij}(d)

    Var[f(x)] ~=       \sum_{ i } \sigma_{ii} f_{,i}^2(d)
               +    2  \sum_{i<j} \sigma_{ij} f_{,i}(d) f_{,j}(d)

which are accurate to second- and first-order, respectively.  Here f_{,i}
denotes partial differentiation with respect to scalar component x_i and
f_{,ij} denotes differentiation with respect to components x_i and x_j.

Taylor Series Methods are discussed at length by Hugh W. Coleman
in "Experimentation, validation, and uncertainty analysis for
engineers." John Wiley & Sons, 3rd edition, 2009. ISBN 0470168889.
The results above are derived in section "Estimating uncertainty in
derived quantities" within Suzerain's perfect gas model document.

Results for Var[f(x)] should be consistent with Table 1 of the article "Notes
on the use of propagation of error formulas." by H. H. Ku appearing in Journal
of Research of the National Bureau of Standards. Section C: Engineering and
Instrumentation, 70C(4):263-273, October 1966.  ISSN 0022-4316.  That article
discusses only the first-order approximation to E[f(x)].
'''
from __future__ import division, print_function
import collections
import fileinput
import itertools
import sympy
import sympy.parsing.sympy_parser
import sys
import tempfile
from sympy.physics.units import Unit

# TODO How to handle uncertainty in derivatives of measured quantities?
# Likely verdict: Chain out derivatives
#                 Deal with correlation between state and derivatives

"Symbolic constants known at parse time to have zero derivatives."
constants = {
    'alpha': Unit('Ratio of bulk to dynamic viscosity', 'alpha'),
    'beta':  Unit('Temperature power law exponent',     'beta' ),
    'gamma': Unit('Ratio of specific heats',            'gamma'),
    'Kn':    Unit('Knudsen number',                     'Kn'   ),
    'Ma':    Unit('Mach number',                        'Ma'   ),
    'Pr':    Unit('Prandtl number',                     'Pr'   ),
    'Re':    Unit('Reynolds number',                    'Re'   ),
}

def parse(f, symbol_table=None):
    r'''
    Given a SymPy expression f or any string parsable as such, produce
    a SymPy expression prepared for further processing by methods
    within this module.  This provides a common extension point for
    injecting known constants (e.g. the Reynolds number Re) and other
    module-specific handling into the parsing process.
    '''
    if isinstance(f, basestring):
        if symbol_table is None:
          t = constants
        else:
          t = constants.copy()
          t.update(symbol_table)
        f = sympy.parsing.sympy_parser.parse_expr(f, t)
    return f

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
                    symbol_table[symbol] = parse(expr, symbol_table)
                except SyntaxError as e:
                    e.filename = fileinput.filename()
                    e.lineno   = fileinput.lineno()
                    raise

    # ...being sure to clean up after our use of fileinput
    finally:
        fileinput.close()

    return symbol_table

def preprocessor(filenames):
    r'''
    Accumulate semicolon-separated statements with C++-style comments.
    TODO Output is a temporary file suitable for input to parser.
    '''
    stmts  = [ [] ]
    # Process input line-by-line...
    try:
        for line in fileinput.input(filenames):

            # ...remove comments defined as the first '//' observed
            head, sep, tail = line.partition('//')
            line = head if head else tail

            # ...and add any statements separated by semis into stmts
            # being careful to permit continuation from prior lines.
            for (i, s) in enumerate(line.split(';')):
                if not i:
                    stmts[-1].append(s.strip())
                else:
                    stmts.append([s.strip()])

    # ...being sure to clean up after our use of fileinput
    finally:
        fileinput.close()

    # Join stmts into list
    return filter(len, [' '.join(s).strip() for s in stmts])


def partials(f):
    r'''
    Given a SymPy expression f or any string parsable as such by parse(),
    produce a defaultdict df where referencing df[x] produces the
    precomputed result f.diff(x).simplify().factor() for any x.

    >>> a, b, c = sympy.symbols('a, b, c')
    >>> df = partials(b**2 + a + 1)
    >>> df[a], df[b], df[c]
    (1, 2*b, 0)

    >>> df = partials("1 + 2 + 3")
    >>> df.keys()
    []

    Fluid-related symbolic constants are treated as genuinely constant:
    >>> df = partials("x + alpha + beta + gamma + Kn + Ma + Pr + Re + y")
    >>> df.keys()
    [x, y]
    '''
    f = parse(f)
    df = collections.defaultdict(lambda: sympy.Integer(0))
    for x in f.free_symbols:
        df[x] = f.diff(x).factor().simplify().factor()
    return df


def mixed_partials(f, df=None):
    r'''
    Given a SymPy expression f or any string parsable as such by parse(),
    produce a defaultdict of defaultdicts ddf where referencing ddf[x][y]
    produces the precomputed result f.diff(x,y).simplify().factor()
    for any x and y.

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
    Given a SymPy expression f or any string parsable as such by parse(),
    return a list wherein unique tuples represents moments necessary
    for computing an estimate of E[f(x)] and Var[f(x)].

    >>> prerequisites('1')
    []

    >>> prerequisites('log(x)')
    [(x,), (x, x)]

    >>> prerequisites('x*y')
    [(x,), (x, x), (x, y), (y,), (y, y)]

    >>> prerequisites('a + x*y')
    [(a,), (a, a), (a, x), (a, y), (x,), (x, x), (x, y), (y,), (y, y)]
    '''
    f = parse(f)

    # Implementation heavily relies on set addition semantics combined
    # with the fact that all derivatives have been precomputed prior
    # to iteration.  Because the caller might know something we do not
    # about he or she wants to compute a subexpression (hinted via
    # df or ddf arguments), we look at all possible terms rather than
    # removing those which can be eliminated by smoothness or symmetry.
    m = set()

    # Quantities necessary to compute first-order Var[f(x)]
    if df is None:
        df = partials(f)
    ## Term:          \sum_{ i } \sigma_{ii} f_{,i}^2(d)
    for i in df.keys():
        if not df[i].is_zero:
            for s in df[i].free_symbols:
                m.add((s,))
            m.add((i, i))                    # Canonical
    ## Term: +    2   \sum_{i<j} \sigma_{ij} f_{,i}(d) f_{,j}(d)
    for (i, j) in itertools.combinations(df.keys(), 2):
        fifj = (df[i] * df[j]).simplify()
        if not fifj.is_zero:
            for s in fifj.free_symbols:
                m.add((s,))
            m.add(tuple(sorted([i, j])))     # Canonicalize

    # Quantities necessary to compute second-order E[f(x)]
    if ddf is None:
        ddf = mixed_partials(f, df)
    ## Term:    f(d)
    for s in f.free_symbols:
        m.add((s,))
    ## Term: + (1/2) \sum_{ i } \sigma_{ii} f_{,ii}(d)
    for i in ddf.keys():
        if not ddf[i][i].is_zero:
            for s in ddf[i][i].free_symbols:
                m.add((s,))
            m.add((i, i))                    # Canonical
    ## Term: +       \sum_{i<j} \sigma_{ij} f_{,ij}(d)
    for (i, j) in itertools.combinations(ddf.keys(), 2):
        if not ddf[i][j].is_zero:
            for s in ddf[i][j].free_symbols:
                m.add((s,))
            m.add(tuple(sorted([i, j])))     # Canonicalize

    return sorted(m)

def expectation(f, ddf=None):
    r'''
    Prepare a map detailing how to compute a second-order approximation
    of E[f(x)].  Keys in the map are either 1 or tuples representing
    covariance scaling factors pre-multiplying the maps' values.
    The maps' values should be evaluated using sample means.

    >>> x = sympy.symbols('x')
    >>> E = expectation(x**2)
    >>> len(E), E[1], E[(x, x)]
    (2, x**2, 1)

    >>> x, y = sympy.symbols('x, y')
    >>> E = expectation(x*y)
    >>> len(E), E[1], E[(x, y)]
    (2, x*y, 1)
    '''
    f = parse(f)
    if ddf is None:
        ddf = mixed_partials(f)

    # Accumulate terms necessary to compute second-order E[f(x)]
    E = collections.defaultdict(lambda: sympy.Integer(0))

    ## Term:    f(d)
    E[1] = f
    ## Term: + (1/2) \sum_{ i } \sigma_{ii} f_{,ii}(d)
    for i in ddf.keys():
        if not ddf[i][i].is_zero:
            E[(i,i)] += (ddf[i][i]/2).simplify()   # Canonical
    ## Term: +       \sum_{i<j} \sigma_{ij} f_{,ij}(d)
    for (i, j) in itertools.combinations(ddf.keys(), 2):
        if not ddf[i][j].is_zero:
            E[tuple(sorted([i, j]))] += ddf[i][j]  # Canonicalize

    return E

def variance(f, df=None):
    r'''
    Prepare a map detailing how to compute a first-order approximation
    of var[f(x)].  Keys in the map are either 1 or tuples representing
    covariances scaling factors pre-multiplying the maps' values.
    The maps' values should be evaluated using sample means.

    >>> x, y = sympy.symbols('x, y')
    >>> Var = variance(2*x + 3*y)
    >>> len(Var), Var[(x, x)], Var[(y, y)], Var[(x, y)]
    (3, 4, 9, 12)

    >>> x, y = sympy.symbols('x, y')
    >>> Var = variance(x/y)
    >>> len(Var), Var[(x, x)], Var[(y, y)], Var[(x, y)]
    (3, y**(-2), x**2/y**4, -2*x/y**3)

    >>> x = sympy.symbols('x')
    >>> Var = variance(x / (1 + x))
    >>> len(Var), Var[(x, x)]
    (1, (x + 1)**(-4))
    '''
    f = parse(f)
    if df is None:
        df = partials(f)

    # Accumulate terms necessary to compute first-order Var[f(x)]
    Var = collections.defaultdict(lambda: sympy.Integer(0))

    ## Term:          \sum_{ i } \sigma_{ii} f_{,i}^2(d)
    for i in df.keys():
        if not df[i].is_zero:
            Var[(i, i)] += (df[i]**2).simplify()   # Canonical
    ## Term: +    2   \sum_{i<j} \sigma_{ij} f_{,i}(d) f_{,j}(d)
    for (i, j) in itertools.combinations(df.keys(), 2):
        twofifj = (2*df[i]*df[j]).simplify()
        if not twofifj.is_zero:
            Var[tuple(sorted([i, j]))] += twofifj  # Canonicalize

    return Var

# def main(args):
#     symbol_table = parser([])
#     return 0
#
# if __name__=='__main__':
#     sys.exit(main(*sys.argv[1:]))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
