#!/usr/bin/env python3
"""Find a docstring that ADVERTISES A CALL THE CODE DOES NOT ACCEPT.  READ-ONLY.

    usage:   python3 tools/checkDocstringKeystrings.py [--all]

WHY THIS EXISTS.  On 07-Nov-2026 six real defects were found during a hygiene run, and EVERY ONE
was found by asking a documentation question rather than by scanning for bugs:

    Basics.LinearMesh()                  returned a GLegenreMesh, silently
    Pulse.GaussianSimplified()           could not construct at all
    ExportAdas.formatExponential         returned 8 characters where its docstring promised 7
    Cascade.simulateRelaxationCurve      claimed an Array where it returns a Tuple
    Basics.read(::ReadMixingFileGrasp18) raised MethodError on every call
    PeriodicTable.getIsotopeData         advertised a ("mass", Z) call that has never existed

FIVE OF THE SIX ARE CALLED BY NOTHING in src/, examples/ or test/, so no test and no coverage
measure could reach them.  What exposed them was having to state in writing what each function
accepts and gives back.  This script mechanises the one part of that which a machine can do:
JAC's string-dispatched functions test their first argument against literal keystrings, and the
docstring lists the keystrings it offers, so the two lists can be compared.

    advertised but NOT TESTED  ->  a documented call that raises.  This is the defect.
    tested but NOT ADVERTISED  ->  an undocumented option.  Reported separately;  often harmless
                                   (a deprecated spelling kept working) and sometimes a real gap.

IT IS DELIBERATELY QUIET ABOUT WHAT IT CANNOT KNOW.  A function whose docstring lists no
keystring is skipped;  so is one whose tests are built from a variable rather than a literal,
since the accepted set is then not in the source at all.  The alternative -- guessing -- is how
every miscount of that same day happened.
"""
import os, re, sys, collections

SRC  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")
SKIP = {"module-RacahAlgebra-inc-sumrules.jl", "module-RacahAlgebra-inc-special.jl"}
OPEN = re.compile(r'^\s*(function|for|while|if|let|begin|try|quote|struct|macro|do)\b')
SHOWALL = "--all" in sys.argv


def blockEnd(L, k):
    depth, i = 0, k
    while i < len(L):
        s = L[i].strip()
        if s.startswith("#") or s.startswith('"""'):  i += 1;  continue
        if OPEN.match(L[i]) and not re.search(r'\bend\b\s*$', L[i]):  depth += 1
        elif re.match(r'^\s*end\b', L[i]):
            depth -= 1
            if depth == 0:  return i
        i += 1
    return min(k+500, len(L)-1)


missing, extra, checked = [], [], 0
family = {}

for fn in sorted(f for f in os.listdir(SRC) if f.endswith(".jl") and f not in SKIP):
    L = open(os.path.join(SRC, fn), errors="replace").read().split("\n")
    i = 0
    while i < len(L):
        if not (L[i].startswith('"""') and L[i].count('"""') == 1):  i += 1;  continue
        j = i+1
        while j < len(L) and '"""' not in L[j]:  j += 1
        k = j+1
        while k < len(L) and L[k].strip() == "":  k += 1
        m = re.match(r'^function\s+([A-Za-z_][\w.!]*)\s*\(([^)]*)', L[k]) if k < len(L) else None
        if not m:  i = j+1;  continue
        # only a function whose FIRST argument is a String is string-dispatched
        first = m.group(2).split(",")[0]
        if "::String" not in first and "::AbstractString" not in first:  i = j+1;  continue
        doc  = "\n".join(L[i:j+1])
        body = "\n".join(L[k:blockEnd(L, k)+1])
        # the keystrings the docstring offers: a quoted literal inside a backticked call, or a
        # quoted literal at the start of a `+ ("...", ...)` continuation line
        adv = set()
        for mm in re.finditer(r'`[^`]*?\(\s*"([^"]+)"', doc):   adv.add(mm.group(1))
        for mm in re.finditer(r'^\s*\+?\s*`?\(\s*"([^"]+)"', doc, re.M):  adv.add(mm.group(1))
        ## THE KEYSTRINGS ARE COLLECTED FROM EVERY METHOD OF THAT NAME IN THE FILE, not from the
        ## one method below the docstring.  JAC writes ONE docstring for a string-dispatched
        ## FAMILY, with `+ ("...", ...)` continuation blocks for the further members, and the
        ## accepted set is then spread over several methods.  Reading only the next body reported
        ## 10 keystrings of `Defaults.setDefaults` as unaccepted when every one of them works.
        tst = set()
        nameOfArg = first.split("::")[0].strip()
        if not nameOfArg:  i = j+1;  continue
        for mk, ml in enumerate(L):
            if not re.match(r'^function\s+' + re.escape(m.group(1)) + r'\s*\(', ml):  continue
            fam = "\n".join(L[mk:blockEnd(L, mk)+1])
            for a in re.finditer(r'\b(\w+)\s*==\s*"([^"]+)"', fam):      tst.add(a.group(2))
            for a in re.finditer(r'"([^"]+)"\s*==\s*\b\w+', fam):        tst.add(a.group(1))
            ## A KEYSTRING MAY ITSELF CONTAIN A BRACKET -- "density: from [g/cm^3] to atomic",
            ## "energy: from wavelength [nm] to atomic" -- so a list body matched as `[^\]]*`
            ## stops at the unit's own closing bracket and loses the entry.  Every quoted string
            ## on a line that opens an `in [` list is taken instead, which is exact for JAC's
            ## one-line `elseif sa in ["a", "b"]` dispatch and reported two accepted keystrings
            ## as unaccepted before it was fixed.
            for ln2 in fam.split("\n"):
                if re.search(r'\b\w+\s+in\s*\[', ln2):
                    tst.update(re.findall(r'"([^"]*)"', ln2))
        if not adv or not tst:  i = j+1;  continue
        ## BOTH SIDES ARE UNIONED OVER THE FAMILY, and that symmetry is the point.  JAC writes one
        ## docstring per string-dispatched FAMILY with `+ ("...", ...)` continuation blocks, and
        ## spreads the accepted keystrings over several methods.  Unioning only the TESTED side
        ## reported every method as carrying dozens of undocumented options -- they were
        ## documented, in a sibling method's block.  So the family is accumulated here and
        ## reported ONCE per function name, below.
        fam = family.setdefault((fn, m.group(1)), {"adv": set(), "tst": set(), "line": k+1})
        fam["adv"].update(adv);  fam["tst"].update(tst)
        fam["line"] = min(fam["line"], k+1)
        i = j+1

for (fn, nm), v in sorted(family.items()):
    checked += 1
    gone = sorted(a for a in v["adv"] if a not in v["tst"])
    more = sorted(t for t in v["tst"] if t not in v["adv"])
    if gone:  missing.append((fn, v["line"], nm, gone))
    if more:  extra.append((fn, v["line"], nm, more))

print("string-dispatched FAMILIES whose docstrings and bodies could both be read: %d" % checked)
print("=" * 96)
print("ADVERTISED BUT NOT ACCEPTED -- a documented call that raises: %d" % len(missing))
for fn, ln, nm, ks in missing:
    print("   %s:%d  %s" % (fn, ln, nm))
    for x in ks:  print('        the docstring offers "%s" and the body never tests for it' % x)
print()
print("ACCEPTED BUT NOT DOCUMENTED -- an undocumented option: %d" % len(extra))
if SHOWALL:
    for fn, ln, nm, ks in extra:
        print("   %s:%d  %s" % (fn, ln, nm))
        print("        %s" % ", ".join('"%s"' % x for x in ks[:8]) + (" ..." if len(ks) > 8 else ""))
else:
    print("   (pass --all to list them;  an undocumented option is often a deprecated spelling)")
