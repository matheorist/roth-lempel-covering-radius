Supplementary exact-computation code

File: supplementary_computations.py

The script uses Python's standard library and exact arithmetic over prime
fields. Running it without command-line options executes four checks tied to
the present manuscript: the degree-gap-threshold MDS example over F_17, the
closed k=3 minimal scan over F_13, the rooted six-set scan over F_13, and the
double-root overlap scan over F_23. The optional
scans independently compare the closed tangent-pair and overlap formulas with
exhaustive finite-field enumeration; all checks terminate with assertions if a
discrepancy is found. For example,

    python supplementary_computations.py --closed-minimal-scan --closed-minimal-q 19
    python supplementary_computations.py --nonzero-branch-example
    python supplementary_computations.py --rooted-six-scan --rooted-six-q 13
    python supplementary_computations.py --double-root-scan --double-root-q 23

Further verification modes are listed by

    python supplementary_computations.py --help

The code is supplementary material and is not required to read the proofs.
