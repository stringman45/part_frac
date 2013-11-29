"""
This module contains various methods with the main purpose to be used by the
part_frac method which returns the partial fraction expansion of a rational
function in one indeterminant.

To run the file use the following command in the SAGE shell.

sage: runfile 'part_frac.sage'

To run the file as a module, you need to preparse the module as a python file;
i.e., use the command in a linux shell (NOT in the SAGE shell):

$ sage --preparse part_frac.sage

This will create a file part_frac.py, which can then be imported as a module
like in normal python.

Next, just call the part_frac method with the desired inputs inside the SAGE
shell. For examples, see the docstrings of the various methods.

Because the purpose of this module is for educational purposes, I used methods
for polynomial division, polynomial GCD, etc., that I wrote instead of Sage's
built-in methods. It should be noted though that I do use Sage's factor
method.

All of the methods in this module are based on exercises in Joel S. Cohen's
book Computer Algebra And Symbolic Computation: Mathematical Methods.

Note that I have not included proper SAGE docstring formatting in this 
module-docstring because I have not found a way to actually see this docstring
in the SAGE command shell.

Name:       John Kluesner
Date:       28 Nov, 2013
Email:      stringman45@gmail.com
"""
def part_frac(f, x):
    r"""
    Return the partial fraction expansion of ``f == u/v``.

    INPUT:

    - ``f`` -- A rational function in the indeterminant ``x``.
    - ``x`` -- A symbol.

    OUTPUT:

    - The partial fraction expansion of ``f == u/v`` in the most general form.

    EXAMPLES:

    The basic functionality is::

        sage: f = (x^3+7*x^2+26*x-105) / (5*x^3-25*x^2-40*x+240)
        sage: w = part_frac(f, x); w
        3/(x - 4) + 5/(x - 4)^2 - 3/5/(x + 3) + 1/5
        sage: bool(w == f)
        True

    Observe that if ``f == u/v`` is not a proper rational function, part_frac
    will divide ``u`` by ``v``. A proper rational function is rational function
    such that the degree of the numerator is less than the degree of the
    denominator::

        sage: u = x^4 + x^3 + x^2 + x + 1
        sage: v = x^2 + 3*x + 2
        sage: u.degree(x) < v.degree(x)
        False
        sage: f = u/v
        sage: w = part_frac(f, x); w
        1/(x + 1) - 11/(x + 2) + x^2 - 2*x + 5
        sage: bool(w == f)
        True

    However, a proper rational function with a denominator that is
    irreducible will just return the rational function::

        sage: u = x^3 + 5*x - 4
        sage: v = x^4 + 1
        sage: f = u/v
        sage: part_frac(f, x)
        (x^3 + 5*x - 4)/(x^4 + 1)

    AUTHORS:

    -John Kluesner (2013-07-12)
    """
    u = f.numerator()
    v = f.denominator()
    if v.degree(x) <= u.degree(x):
        quot, rmdr = poly_div(u, v, x)
        return quot + part_frac_3(rmdr, v.factor(), x)
    else:
        return part_frac_3(u, v.factor(), x)

def poly_div(u,v,x):
    r"""
    Divide the polynomial ``u`` by ``v``.

    INPUT:

    - ``u`` -- A polynomial in ``x``.
    - ``v`` -- A polynomial in ``x``.
    - ``x`` -- A symbol.

    OUTPUT:

    - The quotient and remainder in a list, [quot, rmdr].

    EXAMPLES:

    The basic functionality is::

        sage: u = 5*x^4 + 2*x^3 + 4*x + 1
        sage: v = x^2 + 5*x
        sage: poly_div(u, v, x)
        [5*x^2 - 23*x + 115, -571*x + 1]
        sage: poly_div(x^2 + 5*x + 6, x + 2, x)
        [x + 3, 0]

    AUTHORS:

    -John Kluesner (2013-07-06)
    """
    quot = 0
    rmdr = u
    rmdr_deg = rmdr.degree(x)
    v_deg = v.degree(x)
    lead_coeff_v = v.leading_coeff(x)
    while rmdr_deg >= v_deg and not bool(rmdr == 0):
        lead_coeff_rmdr = rmdr.leading_coeff(x)
        s = lead_coeff_rmdr / lead_coeff_v
        quot = quot + s*x^(rmdr_deg-v_deg)
        rmdr = (rmdr - lead_coeff_rmdr*x^rmdr_deg - (v - lead_coeff_v*x^v_deg)\
                    *s*x^(rmdr_deg-v_deg)).expand()
        rmdr_deg = rmdr.degree(x)
    return [quot, rmdr]


def poly_ext_euc(u, v, x):
    r"""
    Return the greatest common divisor of ``u`` by ``v`` and the polynomials
    ``A`` and ``B`` such that ``Au + Bv = gcd(u ,v)``.

    INPUT:

    - ``u`` -- A polynomial in ``x``.
    - ``v`` -- A polynomial in ``x``.
    - ``x`` -- A symbol.

    OUTPUT:

    - The gcd and the polynomials ``A, B`` as described in the description as
      a list [g, A, B].

    EXAMPLES:

    The basic functionality is::

        sage: u = x^7 - 4*x^5 - x^2 + 4
        sage: v = x^5 - 4*x^3 - x^2 + 4
        sage: g, A, B = poly_ext_euc(u, v, x)

    Here's a test of equality for ``Au + Bv = gcd(u, v)``::

        sage: w = A*u + B*v == g
        sage: bool(w)
        True

    AUTHORS:

    -John Kluesner (2013-07-09)
    """
    if bool(u == 0) and bool(v == 0):
        return [0, 0, 0]
    else:
        poly_1 = u
        poly_2 = v
        cof_poly1_1 = 1
        cof_poly1_2 = 0
        cof_poly2_1 = 0
        cof_poly2_2 = 1
        while not bool(poly_2 == 0):
            quot, rmdr = poly_div(poly_1, poly_2, x)
            cof_poly1 = cof_poly1_1 - quot*cof_poly1_2
            cof_poly2 = cof_poly2_1 - quot*cof_poly2_2
            cof_poly1_1 = cof_poly1_2
            cof_poly1_2 = cof_poly1
            cof_poly2_1 = cof_poly2_2
            cof_poly2_2 = cof_poly2
            poly_1 = poly_2
            poly_2 = rmdr
        lead_coeff = poly_1.leading_coeff(x)
        cof_poly1_1 = (cof_poly1_1/lead_coeff).expand()
        cof_poly2_1 = (cof_poly2_1/lead_coeff).expand()
        poly_1 = (poly_1/lead_coeff).expand()
        return [poly_1, cof_poly1_1, cof_poly2_1]

def poly_expansion(u,v,x,t):
    r"""
    Return the polynomial expansion of ``u`` in ``v``.

    INPUT:

    - ``u`` -- A polynomial in ``x``.
    - ``v`` -- A polynomial in ``x``.
    - ``x`` -- A symbol.
    - ``t`` -- A symbol.

    OUTPUT:

    - The polynomial expansion of ``u`` in ``v`` in terms of ``t``.

    EXAMPLES:

    The basic functionality is::

        sage: u = x^5 + 11*x^4 + 51*x^3 + 124*x^2 + 159*x + 86
        sage: v = x^2 + 4*x + 5
        sage: t = var('t')
        sage: w = poly_expansion(u, v, x, t)
        sage: print(w)
        t^2*x + 3*t^2 + t*x + 2*t + x + 1

    Using a linear polynomial for ``v`` we get all rational coefficients::

        sage: a = 51*x^3 + 124*x^2 + 159*x + 86
        sage: b = 4*x + 5
        sage: c = poly_expansion(a, b, x, t)
        sage: print(c)
        51/64*t^3 - 269/64*t^2 + 1409/64*t - 1191/64

    To see the polynomial expansion in terms of ``v``, use a substitution. 
    Be sure to note automatic simplification::

        sage: c.subs(t = b)
        51/64*(4*x + 5)^3 - 269/64*(4*x + 5)^2 + 1409/16*x + 2927/32

    AUTHORS:

    -John Kluesner (2013-07-06)
    """
    if bool(u == 0):
        return 0
    else:
        quot, rmdr = poly_div(u,v,x)
        return (t*poly_expansion(quot,v,x,t) + rmdr).expand()

def part_frac_1(u, v1, v2, x):
    r"""
    Return the polynomials ``u1, u2`` such that ``u/(v1*v2) = u1/v1 + u2/v2``.

    INPUT:

    - ``u`` -- A polynomial in ``x``.
    - ``v1`` -- A polynomial in ``x``.
    - ``v2`` -- A polynomial in ``x`` with deg(v1*v2,x) > deg(u).
    - ``x`` -- A symbol.

    OUTPUT:

    - ``u1`` and ``u2`` as described above in a list, [u1, u2].

    EXAMPLES:

    The basic functionality is::

        sage: u1, u2 = part_frac_1(3, x, (x+1)^2, x); u1,u2
        (3, -3*x - 6)
        sage: bool(3/x/(x+1)^2 == u1/x + u2/(x+1)^2)
        True

    .. NOTE::

        This is mainly used as a sub-procedure by part_frac.

    AUTHORS:

    -John Kluesner (2013-07-12)
    """
    gcd, A, B = poly_ext_euc(v1, v2, x)
    u2 = poly_div((A*u).expand(), v2, x)[1]
    u1 = poly_div((B*u).expand(), v1, x)[1]
    return [u1, u2]

def part_frac_2(u, v, x):
    r"""
    Return the partial fraction expansion of ``u/v``.

    INPUT:

    - ``u`` -- A polynomial in ``x``.
    - ``v`` -- A polynomial in ``x`` in factored form with deg(u,x) < deg(v,x).
    - ``x`` -- A symbol.

    OUTPUT:

    - The partial fraction expansion of ``u/v``.

    EXAMPLES:

    The basic functionality is::

        sage: w = part_frac_2(x+3, (x^2+3) * (x+1)^2 * (x+5), x); w
        -1/56*(5*x + 3)/(x^2 + 3) + 1/32*(3*x + 7)/(x + 1)^2 - 1/224/(x + 5)
        sage: bool(w == (x+3)/(x^2+3)/(x+5)/(x+1)^2)               
        True

    .. NOTE::

        This is mainly used as a sub-procedure by part_frac. Also, more
        simplification can be done, but that is done in part_fact_3.

    AUTHORS:

    -John Kluesner (2013-07-12)
    """
    if not v.operator() is operator.mul:
        return u/v
    else:
        first_op = v.operands()[0]
        rest_ops = v/first_op
        if not first_op.has(x):
            return 1/first_op*part_frac_2(u, rest_ops, x)
        else:
            u1, w = part_frac_1(u, first_op.expand(), rest_ops.expand(), x)
            return u1/first_op + part_frac_2(w, rest_ops, x)

def part_frac_3(u, v, x):
    r"""
    Return the partial fraction expansion of ``u/v``.

    INPUT:

    - ``u`` -- A polynomial in ``x``.
    - ``v`` -- A polynomial in ``x`` in factored form with deg(u,x) < deg(v,x).
    - ``x`` -- A symbol.

    OUTPUT:

    - The partial fraction expansion of ``u/v`` in the most general form.

    EXAMPLES:

    The basic functionality is::

        sage: w = part_frac_3(12*x^2 + 34*x - 153, 5*(x + 3)*(x - 4)^2, x); w
        3/(x - 4) + 5/(x - 4)^2 - 3/5/(x + 3)
        sage: bool(w == (12*x^2 + 34*x - 153) / 5/(x+3)/(x-4)^2)
        True

    .. NOTE::

        This is mainly used as a sub-procedure by part_frac.

    AUTHORS:

    -John Kluesner (2013-07-12)
    """
    w = part_frac_2(u, v, x)
    if w.operator() is operator.add:
        p = 0
        for op in w.operands():
            if op.denominator().operator() is operator.pow:
                numer = op.numerator().expand()
                base = op.denominator().operands()[0]
                t = var('t')
                expanded_numer = poly_expansion(numer, base, x, t)
                subbed_denom = op.denominator().subs({base: t})
                p0 = (expanded_numer/subbed_denom).expand()
                p += p0.subs(t=base)
            else:
                p += op
        return p
    elif w.operator() is operator.mul and \
                w.denominator().operator() is operator.pow:
        numer = w.numerator().expand()
        base = w.denominator().operands()[0]
        t = var('t')
        expanded_numer = poly_expansion(numer, base, x, t).subs(t=base)
        return (expanded_numer/w.denominator()).expand()
    else:
        return w
