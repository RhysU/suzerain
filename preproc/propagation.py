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

which are accurate to second- and first order, respectively.  Here f_{,i}
denotes partial differentiation with respect to scalar component x_i and
f_{,ij} denotes differentiation with respect to components x_i and x_j.
Notice that the first order approximation to E[f(x)] is nothing but f(d).

Taylor Series Methods are discussed at length by Hugh W. Coleman
in "Experimentation, validation, and uncertainty analysis for
engineers." John Wiley & Sons, 3rd edition, 2009. ISBN 0470168889.
The results above are derived in section "Estimating uncertainty in
derived quantities" within Suzerain's perfect gas model document.
'''
from __future__ import division, print_function
import collections
import distutils.version
import fileinput
import itertools
import sympy
import sympy.parsing.sympy_parser
import sys
from sympy.core.function import AppliedUndef

# FIXME Hand check correct behavior against several known Coleman cases


# Separator between a quantity name and the partial derivatives requested. For
# example, when 'foo' is some quantity then 'foo_xyz' is its partial derivative
# with respect to x, y, and z.
daffsep = '__'


def daff(expr, *symbols, **kwargs):
    r'''
    Write Derivative(f(x), x) to be f__x(x), effectively flattening
    it. Named as a contraction of 'derivative affix', the functionality
    is meant to facilitate a derivative naming convention like foo__y(y).

    First and second derivatives are named by per the convention:

    >>> f, x, y = sympy.Function('f'), sympy.Symbol('x'), sympy.Symbol('y')
    >>> daff( (f(x)).diff(x) )
    f__x(x)
    >>> daff( (f(x)).diff(x).diff(x) )
    f__xx(x)

    Regular function application is unaffected:

    >>> daff( sympy.cos(3*f(x, y)) )
    cos(3*f(x, y))

    Multivariate functions behavior assumes partial differentiation commutes:

    >>> daff( (f(x,y)).diff(x).diff(y) )
    f__xy(x, y)
    >>> daff( (f(y,x)).diff(y).diff(x) )
    f__xy(y, x)

    Because calling sympy.diff then daff is cumbersome, a multiargument
    invocation automatically wraps daff around sympy.diff:

    >>> daff(f(x), x)
    f__x(x)
    >>> daff(f(y,x), y, x)
    f__xy(y, x)
    >>> daff(f(x), x, 5)
    f__xxxxx(x)

    '''

    def helper(f, *wrt):
        head, sep, tail = type(f).__name__.partition(daffsep)
        deriv = list(tail)
        deriv.extend(sym.name for sym in wrt)
        name = [head, daffsep]
        name.extend(sorted(deriv))
        return sympy.Function(''.join(name))(*list(f.args))

    if symbols:
        expr = sympy.diff(expr, *symbols, **kwargs).simplify()

    return expr.replace(sympy.Derivative, helper)


def parse(expr, symbol_table=None):
    r'''
    Given a SymPy expression expr or any string parsable as such, produce a
    SymPy expression.  This method provides a common extension point for
    injecting module-specific handling into the parsing process.
    '''
    if isinstance(expr, basestring):
        expr = sympy.parsing.sympy_parser.parse_expr(expr, symbol_table)

    return expr


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


class symboltranscript(collections.OrderedDict):
    r'''
    A collections.OrderedDict subclass for symbol -> SymPy expression
    entries.  Accordingly, entry insertion order is preserved.  Lookup of
    missing entries like 'foo__y' will automatically invoke 'daff(foo,y)',
    insert the result, and return the derivative.
    '''
    def __missing__(self, key):
        'Uses daff(...) to generate missing keys whenever possible.'
        head, sep, tail = key.partition(daffsep)
        if sep and tail:
            if head in self:
                val = daff(self[head], *list(tail))  # Generate derivative...
                self[key] = val                      # ...and insert entry
                return val
            else:
                raise KeyError('%s depends on missing key %s' % (key, head))
        else:
            raise KeyError(key)

    def __contains__(self, key):
        'Augments containment to include keys computable via daff(...).'
        if super(symboltranscript, self).__contains__(key):
            return True
        else:
            head, sep, tail = key.partition(daffsep)
            if sep and tail:
                return super(symboltranscript, self).__contains__(head)
            else:
                return False


def canonical(*exprs):
    r'''
    Produce a canonically ordered tuple of the provided SymPy expressions.
    Older SymPy permitted tuple(sorted([f(x),g(x)])) but that breaks on 0.7.4
    (refer to http://stackoverflow.com/questions/24093363/ for more details).
    Workaround employs SymPy compare member method via sorted's key argument.

    Canonical sorting of symbols:

    >>> x, y = map(sympy.Symbol, 'xy')
    >>> canonical(x, y)
    (x, y)
    >>> canonical(y, x)
    (x, y)
    >>> canonical(x, x)
    (x, x)

    Canonical sorting of functions of the same argument:

    >>> f, g = map(sympy.Function, 'fg')
    >>> canonical(f(x), g(x))
    (f(x), g(x))
    >>> canonical(g(x), f(x))
    (f(x), g(x))

    Canonical sorting of one function with different arguments:

    >>> canonical(f(x), f(y))
    (f(x), f(y))

    >>> canonical(f(y), f(x))
    (f(x), f(y))

    >>> canonical(f(x), f(x))
    (f(x), f(x))

    Canonical sorting of one function invoked with differing arity:

    >>> canonical(f(x, y), f(x))
    (f(x), f(x, y))

    '''
    # Similar to functools.cmp_to_key() but without intermediate cmp function
    class Key(object):
        def __init__(self, obj, *args):
            self.obj = obj

        def __lt__(self, other):
            return self.obj.compare(other.obj) < 0

        def __gt__(self, other):
            return self.obj.compare(other.obj) > 0

        def __eq__(self, other):
            return self.obj.compare(other.obj) == 0

        def __le__(self, other):
            return self.obj.compare(other.obj) <= 0

        def __ge__(self, other):
            return self.obj.compare(other.obj) >= 0

        def __ne__(self, other):
            return self.obj.compare(other.obj) != 0
    return tuple(sorted(exprs, key=Key))


def parser(statement_tuples):
    r'''
    Parse statements from (filename, lineno, statement) tuples into a
    symboltranscript.  Either statements_by_newline() or
    statements_by_semicolon() may be used to generate tuples from input files.

    >>> parser([ ("test", 1, "a=1"     )    # Simple assignment
    ...        , ("test", 2, "b  = a+1")    # Reuse earlier definition
    ...        , ("test", 3, "c =  d+e") ]) # Purely symbolic result okay
    symboltranscript([('a', 1), ('b', 2), ('c', d + e)])

    Lacking assignment, a target name will be generated from the line number:

    >>> parser([ ("test", 4, "f") ])
    symboltranscript([('line4', f)])

    Symbol re-definitions cause LookupErrors identifying the location:

    >>> parser([ ('somefile', 1, 'a=b')
    ...        , ('somefile', 2, 'a=c') ])
    Traceback (most recent call last):
        ...
    LookupError: Symbol 'a' redefined at somefile:2
    '''
    # Accumulate symbol definitions maintaining declaration order
    symbol_table = symboltranscript()

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
    produce a defaultdict df where referencing df[x] produces the result
    f.diff(fx).simplify().factor() for every fx in f.atoms(AppliedUndef).


    >>> f, g, h = map(sympy.Function, 'fgh')
    >>> x, y, z = map(sympy.Symbol,   'xyz')
    >>> log = sympy.log
    >>> df = partials(g(x)**2 + f(x) + 1 + log(h(x)))
    >>> df[f(x)], df[g(x)], df[h(x)], df[log(h(x))]
    (1, 2*g(x), 1/h(x), 0)

    Non-functions are not considered for the list of partial derivatives:

    >>> df = partials("x + y + z + 1 + 2 + 3")
    >>> df.keys()
    []
    '''
    f = parse(f)
    df = collections.defaultdict(lambda: sympy.Integer(0))
    for fx in f.atoms(AppliedUndef):
        df[fx] = f.diff(fx).factor().simplify().factor()
    return df


def mixed_partials(f, df=None):
    r'''
    Given a SymPy expression f or any string parsable as such by parse(),
    produce a defaultdict of defaultdicts ddf where referencing ddf[fx][fy]
    produces the precomputed result f.diff(fx,fy).simplify().factor()
    for all fx and fy in expr.atoms(AppliedUndef).

    >>> f, g, h = map(sympy.Function, 'fgh')
    >>> x, y, z = map(sympy.Symbol,   'xyz')
    >>> ddf = mixed_partials(f(x)**2 + f(x)*g(x) + g(x)**2 + 1)
    >>> ddf[f(x)][f(x)], ddf[f(x)][g(x)], ddf[g(x)][f(x)], ddf[g(x)][g(x)]
    (2, 1, 1, 2)
    >>> ddf[f(x)][h(x)], ddf[h(x)][f(x)], ddf[h(x)][h(x)]
    (0, 0, 0)
    '''
    ddf = collections.defaultdict(
        lambda: collections.defaultdict(lambda: sympy.Integer(0)))
    if df is None:
        df = partials(f)
    for (x, dfdx) in df.iteritems():
        ddf[x] = partials(dfdx)

    return ddf


def prerequisites(f, df=None, ddf=None, order=2):
    r'''
    Given a SymPy expression f or any string parsable as such by parse(),
    return a set wherein unique tuples represents moments necessary to
    estimate E[f(x)] to first- or second order and Var[f(x)] to first order.

    >>> sorted(prerequisites('1 + x*y + log(x/y)'))
    []

    >>> sorted(prerequisites('f(x)*g(y) + a'))
    [(f(x),), (f(x), f(x)), (f(x), g(y)), (g(y),), (g(y), g(y))]

    >>> sorted(prerequisites('f(x)*g(y) + a', order=2))
    [(f(x),), (f(x), f(x)), (f(x), g(y)), (g(y),), (g(y), g(y))]
    '''
    f = parse(f)

    # Implementation heavily relies on set addition semantics combined
    # with the fact that all derivatives have been precomputed prior
    # to iteration.  Because the caller might know something we do not
    # about how he or she wants to compute a subexpression (hinted via
    # df or ddf arguments), we look at all possible terms rather than
    # removing those which can be eliminated by smoothness or symmetry.
    m = set()

    # Quantities necessary to compute first order Var[f(x)]
    if df is None:
        df = partials(f)
    # Term:          \sum_{ i } \sigma_{ii} f_{,i}^2(d)
    for i in df.keys():
        if not df[i].is_zero:
            for s in df[i].atoms(AppliedUndef):
                m.add((s,))
            m.add((i, i))                    # Canonical
    # Term: +    2   \sum_{i<j} \sigma_{ij} f_{,i}(d) f_{,j}(d)
    for (i, j) in itertools.combinations(df.keys(), 2):
        fifj = (df[i] * df[j]).simplify()
        if not fifj.is_zero:
            for s in fifj.atoms(AppliedUndef):
                m.add((s,))
            m.add(tuple(sorted([i, j])))     # Canonicalize

    # Quantities necessary to compute first- or second order E[f(x)]
    assert order == 1 or order == 2
    if ddf is None:
        ddf = mixed_partials(f, df)
    # Term:    f(d)
    for s in f.atoms(AppliedUndef):
        m.add((s,))
    if order > 1:
        # Term: + (1/2) \sum_{ i } \sigma_{ii} f_{,ii}(d)
        for i in ddf.keys():
            if not ddf[i][i].is_zero:
                for s in ddf[i][i].atoms(AppliedUndef):
                    m.add((s,))
                m.add((i, i))                    # Canonical
        # Term: +       \sum_{i<j} \sigma_{ij} f_{,ij}(d)
        for (i, j) in itertools.combinations(ddf.keys(), 2):
            if not ddf[i][j].is_zero:
                for s in ddf[i][j].atoms(AppliedUndef):
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


def expectation(f, ddf=None, order=2):
    r'''
    Prepare momentdict detailing an approximation to E[f(x)].

    Second order approximations are computed by default:

    >>> f, x = sympy.Function('f'), sympy.Symbol('x')
    >>> E = expectation(f(x)**2)
    >>> len(E), E[1], E[(f(x), f(x))]
    (2, f(x)**2, 1)

    >>> g = sympy.Function('g')
    >>> E = expectation(f(x)*g(x))
    >>> len(E), E[1], E[(f(x), g(x))]
    (2, f(x)*g(x), 1)

    First order approximations may be explicitly requested:

    >>> g = sympy.Function('g')
    >>> E = expectation(f(x)*g(x), order=1)
    >>> len(E), E[1], E[(f(x), g(x))]
    (1, f(x)*g(x), 0)

    Non-functions are assumed to lack interesting correlations:

    >>> E = expectation(x**2 * f(x))
    >>> len(E), E[1]
    (1, x**2*f(x))
    '''
    f = parse(f)
    if ddf is None:
        ddf = mixed_partials(f)

    # Accumulate terms necessary to compute first- or second order E[f(x)]
    assert order == 1 or order == 2
    E = momentdict()

    # Term:    f(d)
    E[1] = f
    if order > 1:
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

    Results for Var[f(x)] should be consistent with Table 1 of the
    article "Notes on the use of propagation of error formulas." by
    H. H. Ku appearing in Journal of Research of the National Bureau of
    Standards. Section C: Engineering and Instrumentation, 70C(4):263-273,
    October 1966.  ISSN 0022-4316.  Beware, however, that article discusses
    only the first-order approximation to E[f(x)].

    >>> f, g, x = sympy.Function('f'), sympy.Function('g'), sympy.Symbol('x')
    >>> Var = variance(2*f(x) + 3*g(x))
    >>> len(Var), Var[(f(x), f(x))], Var[(g(x), g(x))], Var[(f(x), g(x))]
    (3, 4, 9, 12)

    >>> Var = variance(f(x)/g(x))
    >>> len(Var), Var[(f(x), f(x))], Var[(g(x), g(x))], Var[(f(x), g(x))]
    (3, g(x)**(-2), f(x)**2/g(x)**4, -2*f(x)/g(x)**3)

    >>> Var = variance(f(x) / (1 + f(x)))
    >>> len(Var), Var[(f(x), f(x))]
    (1, (f(x) + 1)**(-4))

    Non-functions are assumed to lack interesting correlations:

    >>> Var = variance(x**2 * f(x))
    >>> len(Var), Var[(f(x), f(x))]
    (1, x**4)
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
    p.add_argument('-d', '--decl', action='append', default=[],
                   help='file containing SymPy-based declarations defaulting'
                        ' to standard input when declarations required')

    # Control C vs Fortran vs SymPy (default) output for expressions
    g = p.add_mutually_exclusive_group()
    g.add_argument('-c', '--ccode', dest='style', default=str,
                   help='output symbolic expressions as C code',
                   action='store_const',
                   const=lambda expr: sympy.ccode(expr, human=False)[2])
    g.add_argument('-f', '--fcode', dest='style', default=SUPPRESS,
                   help='output symbolic expressions as Fortran code',
                   action='store_const',
                   const=lambda expr: sympy.fcode(expr,
                                                  source_format='free',
                                                  human=False)[2])

    # Control order of the approximation to any used expectations
    g = p.add_mutually_exclusive_group()
    g.add_argument('-1', dest='order', action='store_const', const=1,
                    help='employ 1st order approximation to the expectation')
    g.add_argument('-2', dest='order', action='store_const', const=2,
                    help='employ 2nd order approximation to the expectation')
    p.set_defaults(order=2)

    # Add command-specific subparsers
    sp = p.add_subparsers(title='subcommands, each possibly using declarations',
                          metavar='')
    sp_chk = sp.add_parser('chk', help='run verification sanity checks',
                           description='run verification sanity checks')
    sp_dec = sp.add_parser('dec', help='output known declarations',
                           description='output known declarations')
    sp_lib = sp.add_parser('lib', help='list undefined symbols in declarations',
                           description=prerequisites.__doc__,
                           formatter_class=RawDescriptionHelpFormatter)
    sp_pre = sp.add_parser('pre', help='list prerequisites for E[f(x)], Var[f(x)]',
                           description=prerequisites.__doc__,
                           formatter_class=RawDescriptionHelpFormatter)
    sp_exp = sp.add_parser('exp', help='tabulate terms in E[f(x)]',
                           description=expectation.__doc__,
                           formatter_class=RawDescriptionHelpFormatter)
    sp_var = sp.add_parser('var', help='tabulate terms in Var[f(x)]',
                           description=variance.__doc__,
                           formatter_class=RawDescriptionHelpFormatter)

    # Some of the subparsers take one or more symbols to process
    f_help = 'quantities of interest; if empty, process all declarations'
    sp_dec.add_argument('f', nargs='*', help=f_help)
    sp_lib.add_argument('f', nargs='*', help=f_help)
    sp_pre.add_argument('f', nargs='*', help=f_help)
    sp_exp.add_argument('f', nargs='*', help=f_help)
    sp_var.add_argument('f', nargs='*', help=f_help)

    # Each command dispatches to one of the following methods
    sp_chk.set_defaults(command=command_chk)
    sp_dec.set_defaults(command=command_dec)
    sp_lib.set_defaults(command=command_lib)
    sp_pre.set_defaults(command=command_pre)
    sp_exp.set_defaults(command=command_exp)
    sp_var.set_defaults(command=command_var)

    # Parse the incoming command line
    args = p.parse_args()

    # Supply "all known declarations" behavior for any f arguments first
    # parsing any requested files into one unified dictionary "args.syms".
    # Additionally, provide nicer error messages on unknown symbols
    # (otherwise a messy stacktrace appears and folks doubt the code).
    if ('f' in args):
        args.syms = parser(args.statements(args.decl))
        if not args.f:
            args.f = args.syms.keys()
        else:
            unknown = filter(lambda k: k not in args.syms, args.f)
            if unknown:
                raise LookupError("Requested but not in declarations: %s"
                                  % ', '.join(unknown))

    # Dispatch to the chosen command
    retval = args.command(args)

    # Warn if the results are likely bogus due to an older SymPy version
    if (   distutils.version.StrictVersion(sympy.__version__)
         < distutils.version.StrictVersion('0.7.2')           ):
        print('WARN: %s notes SymPy %s older than minimum 0.7.2'
              % (argv[0], sympy.__version__),
              file=sys.stderr)

    return retval


def command_chk(args):
    r'''Process the 'chk' command on behalf of main()'''
    from doctest import testmod
    failure_count, test_count = testmod(verbose=args.verbosity)
    return failure_count


def command_dec(args):
    r'''Process the 'dec' command on behalf of main()'''
    for qoi in args.f:
        print(qoi, '=', args.style(args.syms[qoi]))
    return 0


def command_lib(args):
    r'''Process the 'lib' command on behalf of main()'''
    freesyms = set()
    freefunc = set()
    for qoi in args.f:
        decl = args.syms[qoi]
        freefunc.update(decl.atoms(AppliedUndef))
        freesyms.update(decl.free_symbols)
    map(print, sorted(freefunc))
    print()
    map(print, sorted(freesyms))


def command_pre(args):
    r'''Process the 'pre' command on behalf of main()'''
    prereqs = set()
    for qoi in args.f:
        prereqs.update(prerequisites(args.syms[qoi], order=args.order))
    for prereq in sorted(prereqs):
        print('E', list(prereq), sep='')


def command_exp(args):
    r'''Process the 'exp' command on behalf of main()'''
    for qoi in args.f:
        m = expectation(args.syms[qoi], order=args.order)
        print('E[', qoi, '] = (\n', m.__str__(args.style), ')\n', sep='')


def command_var(args):
    r'''Process the 'var' command on behalf of main()'''
    for qoi in args.f:
        m = variance(args.syms[qoi])
        print('Var[', qoi, '] = (\n', m.__str__(args.style), ')\n', sep='')


if __name__ == "__main__":
    sys.exit(main(sys.argv))
