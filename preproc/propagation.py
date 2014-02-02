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
Assume also \vec{\epsilon} is a zero-mean, normally-distributed random
variable with some known covariance matrix \Sigma containing scalar
components \sigma_{ij}.

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
Instrumentation, 70C(4):263-273, October 1966.  ISSN 0022-4316.  Beware,
however, that article discusses only the first-order approximation to E[f(x)].
'''
from __future__ import division, print_function
import collections
import fileinput
import itertools
import sympy
import sympy.parsing.sympy_parser
from sympy.physics.units import Unit

# TODO How to handle uncertainty in derivatives of measured quantities?
# Likely verdict: Chain out derivatives
#                 Deal with correlation between state and derivatives

"Symbolic constants known at parse time to have zero derivatives."
constants = {
    'alpha': Unit('Ratio of bulk to dynamic viscosity', 'alpha'),
    'beta':  Unit('Temperature power law exponent',     'beta'),
    'gamma': Unit('Ratio of specific heats',            'gamma'),
    'Kn':    Unit('Knudsen number',                     'Kn'),
    'Ma':    Unit('Mach number',                        'Ma'),
    'Pr':    Unit('Prandtl number',                     'Pr'),
    'Re':    Unit('Reynolds number',                    'Re'),
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


# TODO Line continuation via trailing backslash
def statements_by_newline(files=None):
    r'''
    Generate (filename, lineno, statement) tuples by parsing the provided
    filenames with newline-separated, whitespace-trimmed statements.
    Comments are introduced by a '#' and extend until the end of line.

    >>> from tempfile import NamedTemporaryFile
    >>> with NamedTemporaryFile() as f:
    ...     print("""a=1       # Trailing comments
    ...                        # Not every line must have a statement
    ...              f         # Nor every line involve assignment
    ...           """, file=f)
    ...     f.flush()
    ...     for (_, lineno, stmt) in statements_by_newline(f.name):
    ...         print(lineno, stmt)
    1 a=1
    3 f
    '''
    # Process input line-by-line...
    f = fileinput.FileInput(files)
    for line in f:

        # ...remove comments occurring after the first '#' character
        line, _, _ = line.partition('#')

        # ...trim then yield statement only on nontrivial line
        line = line.strip()
        if line:
            yield (f.filename(), f.filelineno(), line)


# TODO Behavior on lingering statement content without semicolon
def statements_by_semicolon(files=None):
    r'''
    Generate (filename, lineno, statement) tuples by parsing the provided
    filenames with semicolon-separated, whitespace-trimmed statements.
    Comments are introduced by a '//' and extend until the end of line.

    >>> from tempfile import NamedTemporaryFile
    >>> with NamedTemporaryFile() as f:
    ...     print("""a=1;      // Trailing comments may include ';'
    ...              b =       // Statements may span lines
    ...                  c;
    ...              1;2;;     // Multiple may appear with empty ignored
    ...           """, file=f)
    ...     f.flush()
    ...     for (_, lineno, stmt) in statements_by_semicolon(f.name):
    ...         print(lineno, stmt)
    1 a=1
    3 b = c
    4 1
    4 2
    '''
    # Process input line-by-line maintaining any active statement...
    f = fileinput.FileInput(files)
    stmt = []
    for line in f:

        # ...remove comments defined as the first '//' observed
        line, _, _ = line.partition('//')

        # ...and yield any statements separated by semicolons
        # being careful to permit continuation from prior lines.
        while line:
            head, sep, line = line.partition(';')
            head = head.strip()
            if head:
                stmt.append(head)
            if sep and stmt:
                yield (f.filename(), f.filelineno(), ' '.join(stmt))
                del stmt[:]


def parser(statement_tuples):
    r'''
    Parse statements from (filename, lineno, statement) tuples into
    a collections.OrderedDict of symbol -> SymPy expression entries.
    Either statements_by_newline() or statements_by_semicolon() may be used to
    generate tuples from input files.

    >>> parser([ ("test", 1, "a=1"     )    # Simple assignment
    ...        , ("test", 2, "b  = a+1")    # Reuse earlier definition
    ...        , ("test", 3, "c =  d+e") ]) # Purely symbolic result okay
    OrderedDict([('a', 1), ('b', 2), ('c', d + e)])

    Lacking assignment, a target name will be generated from the line number:

    >>> parser([ ("test", 4, "f") ])
    OrderedDict([('line4', f)])

    Symbol re-definitions cause LookupErrors identifying the location:

    >>> parser([ ('somefile', 1, 'a=b')
    ...        , ('somefile', 2, 'a=c') ])
    Traceback (most recent call last):
        ...
    LookupError: Symbol 'a' redefined at somefile:2
    '''
    # Accumulate symbol definitions maintaining declaration order
    symbol_table = collections.OrderedDict()

    # Process each trimmed, comment-less statement...
    for (filename, lineno, stmt) in statement_tuples:

        # ...split lines into "token = rhs" or just "rhs" with the
        # latter implicitly naming the result like "line1234".
        symbol, sep, expr = stmt.partition('=')
        expr = expr.strip()
        if not expr:
            symbol, expr = 'line' + repr(lineno), symbol
        symbol = symbol.strip()
        if symbol in symbol_table:
            raise LookupError("Symbol '%s' redefined at %s:%d"
                              % (symbol, filename, lineno))

        # ...parse and save the result into the known symbol table
        # augmenting any parsing errors with location information
        try:
            symbol_table[symbol] = parse(expr, symbol_table)
        except SyntaxError as e:
            e.filename = filename
            e.lineno = lineno
            raise

    return symbol_table


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
        lambda: collections.defaultdict(lambda: sympy.Integer(0)))
    if df is None:
        df = partials(f)
    for (x, dfdx) in df.iteritems():
        ddf[x] = partials(dfdx)

    return ddf


def prerequisites(f, df=None, ddf=None):
    r'''
    Given a SymPy expression f or any string parsable as such by parse(),
    return a set wherein unique tuples represents moments necessary
    for computing an estimate of E[f(x)] and Var[f(x)].

    >>> sorted(prerequisites('1'))
    []

    >>> sorted(prerequisites('log(x)'))
    [(x,), (x, x)]

    >>> sorted(prerequisites('x*y'))
    [(x,), (x, x), (x, y), (y,), (y, y)]

    >>> sorted(prerequisites('a + x*y'))
    [(a,), (a, a), (a, x), (a, y), (x,), (x, x), (x, y), (y,), (y, y)]
    '''
    f = parse(f)

    # Implementation heavily relies on set addition semantics combined
    # with the fact that all derivatives have been precomputed prior
    # to iteration.  Because the caller might know something we do not
    # about how he or she wants to compute a subexpression (hinted via
    # df or ddf arguments), we look at all possible terms rather than
    # removing those which can be eliminated by smoothness or symmetry.
    m = set()

    # Quantities necessary to compute first-order Var[f(x)]
    if df is None:
        df = partials(f)
    # Term:          \sum_{ i } \sigma_{ii} f_{,i}^2(d)
    for i in df.keys():
        if not df[i].is_zero:
            for s in df[i].free_symbols:
                m.add((s,))
            m.add((i, i))                    # Canonical
    # Term: +    2   \sum_{i<j} \sigma_{ij} f_{,i}(d) f_{,j}(d)
    for (i, j) in itertools.combinations(df.keys(), 2):
        fifj = (df[i] * df[j]).simplify()
        if not fifj.is_zero:
            for s in fifj.free_symbols:
                m.add((s,))
            m.add(tuple(sorted([i, j])))     # Canonicalize

    # Quantities necessary to compute second-order E[f(x)]
    if ddf is None:
        ddf = mixed_partials(f, df)
    # Term:    f(d)
    for s in f.free_symbols:
        m.add((s,))
    # Term: + (1/2) \sum_{ i } \sigma_{ii} f_{,ii}(d)
    for i in ddf.keys():
        if not ddf[i][i].is_zero:
            for s in ddf[i][i].free_symbols:
                m.add((s,))
            m.add((i, i))                    # Canonical
    # Term: +       \sum_{i<j} \sigma_{ij} f_{,ij}(d)
    for (i, j) in itertools.combinations(ddf.keys(), 2):
        if not ddf[i][j].is_zero:
            for s in ddf[i][j].free_symbols:
                m.add((s,))
            m.add(tuple(sorted([i, j])))     # Canonicalize

    return m


class momentdict(collections.defaultdict):
    r'''
    A map detailing how to compute some approximate moment like E[f(x)]
    or Var[f(x)].  Keys in the map are either 1 or tuples representing
    covariance scaling factors pre-multiplying the maps' values.
    The maps' values should be evaluated using sample means.
    '''
    def __init__(self):
        super(momentdict, self).__init__(lambda: sympy.Integer(0))

    def __str__(self, sym=str):
        r'''
        Produce string representation using sym to format Sympy expressions
        '''
        s = []
        p = '  '
        if 1 in self:
            if len(self) > 1:
                s.append('  ')
                p = '+ '
            s.append('(')
            s.append(sym(self[1]))
            s.append(')\n')
        for term in sorted(k for k in self.keys() if k != 1):
            s.append(p)
            s.append('E')
            s.append(str(list(term)))
            s.append(' * (')
            s.append(sym(self[term]))
            s.append(')\n')
            p = '+ '
        return ''.join(s)

def expectation(f, ddf=None):
    r'''
    Prepare momentdict detailing the second-order approximation to E[f(x)].

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
    E = momentdict()

    # Term:    f(d)
    E[1] = f
    # Term: + (1/2) \sum_{ i } \sigma_{ii} f_{,ii}(d)
    for i in ddf.keys():
        if not ddf[i][i].is_zero:
            E[(i, i)] += (ddf[i][i] / 2).simplify()  # Canonical
    # Term: +       \sum_{i<j} \sigma_{ij} f_{,ij}(d)
    for (i, j) in itertools.combinations(ddf.keys(), 2):
        if not ddf[i][j].is_zero:
            E[tuple(sorted([i, j]))] += ddf[i][j]    # Canonicalize

    return E


def variance(f, df=None):
    r'''
    Prepare momentdict detailing the first-order approximation to Var[f(x)].

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
    Var = momentdict()

    # Term:          \sum_{ i } \sigma_{ii} f_{,i}^2(d)
    for i in df.keys():
        if not df[i].is_zero:
            Var[(i, i)] += (df[i] ** 2).simplify()  # Canonical
    # Term: +    2   \sum_{i<j} \sigma_{ij} f_{,i}(d) f_{,j}(d)
    for (i, j) in itertools.combinations(df.keys(), 2):
        twofifj = (2 * df[i] * df[j]).simplify()
        if not twofifj.is_zero:
            Var[tuple(sorted([i, j]))] += twofifj   # Canonicalize

    return Var


def main(argv):
    # Define arguments applicable to all commands
    from argparse import ArgumentParser, RawDescriptionHelpFormatter, SUPPRESS
    p = ArgumentParser(description=__doc__,
                       formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('-v', '--verbosity', action='count',
                   help='increase verbosity; may be supplied repeatedly')
    p.add_argument('-n', '--newline', action='store_const', dest='statements',
                   const=statements_by_newline,
                   default=statements_by_semicolon,
                   help='expect newline-delimited input with "#" comments;'
                        ' otherwise semicolons delimit with "//" comments')
    p.add_argument('-d', '--decl', type=str,
                   help='file containing SymPy-based declarations;'
                        ' process standard input when "-" supplied')
    g = p.add_mutually_exclusive_group()
    g.add_argument('-c', '--ccode', dest='repr', default=str,
                   help='output symbolic expressions as C code',
                   action='store_const', const=sympy.ccode)
    g.add_argument('-f', '--fcode', dest='repr', default=SUPPRESS,
                   help='output symbolic expressions as Fortran code',
                   action='store_const', const=sympy.fcode)

    # Add command-specific subparsers
    sp = p.add_subparsers(title='Operations to perform on declarations',
                          help='Exactly one operation must be supplied')
    sp_chk = sp.add_parser('chk', help='Run verification sanity checks',
                           description='Run verification sanity checks')
    sp_dec = sp.add_parser('dec', help='Output known declarations',
                           description='Output known declarations')
    sp_pre = sp.add_parser('pre', help='Prerequisites for E[f(x)], Var[f(x)]',
                           description=prerequisites.__doc__,
                           formatter_class=RawDescriptionHelpFormatter)
    sp_exp = sp.add_parser('exp', help='Tabulate terms in E[f(x)]',
                           description=expectation.__doc__,
                           formatter_class=RawDescriptionHelpFormatter)
    sp_var = sp.add_parser('var', help='Tabulate terms in Var[f(x)]',
                           description=variance.__doc__,
                           formatter_class=RawDescriptionHelpFormatter)

    # Some of the subparsers take one or more symbols to process
    f_help='quantities of interest; if empty, process all declarations'
    sp_dec.add_argument('f', nargs='*', help=f_help)
    sp_pre.add_argument('f', nargs='*', help=f_help)
    sp_exp.add_argument('f', nargs='*', help=f_help)
    sp_var.add_argument('f', nargs='*', help=f_help)

    # Each command dispatches to one of the following methods
    sp_chk.set_defaults(command=command_chk)
    sp_dec.set_defaults(command=command_dec)
    sp_pre.set_defaults(command=command_pre)
    sp_exp.set_defaults(command=command_exp)
    sp_var.set_defaults(command=command_var)

    # Parse the incoming command line
    args = p.parse_args()

    # Parse any requested files into one unified symbol dictionary
    syms = parser(args.statements(args.decl) if args.decl else [])

    # Supply "all known declarations" behavior for any f arguments
    if ('f' in args) and (not args.f):
        args.f = syms.keys()

    # Dispatch to the chosen command passing arguments and symbols
    return args.command(args, syms)


def command_chk(args, syms):
    r'''Process the 'chk' command on behalf of main()'''
    from doctest import testmod
    failure_count, test_count = testmod(verbose=args.verbosity)
    return failure_count


def command_dec(args, syms):
    r'''Process the 'dec' command on behalf of main()'''
    for qoi in args.f:
        print(qoi, '=', syms[qoi])
    return 0


def command_pre(args, syms):
    r'''Process the 'pre' command on behalf of main()'''
    prereqs = set()
    for qoi in args.f:
        prereqs.update(prerequisites(syms[qoi]))
    for prereq in sorted(prereqs):
        print('E', list(prereq), sep='')


def command_exp(args, syms):
    r'''Process the 'exp' command on behalf of main()'''
    for qoi in args.f:
        m = expectation(syms[qoi])
        print('E[', qoi, '] = (\n', m.__str__(args.repr), ')\n', sep='')


def command_var(args, syms):
    r'''Process the 'var' command on behalf of main()'''
    for qoi in args.f:
        m = variance(syms[qoi])
        print('Var[', qoi, '] = (\n', m.__str__(args.repr), ')\n', sep='')


if __name__ == "__main__":
    from sys import argv, exit
    exit(main(argv))
