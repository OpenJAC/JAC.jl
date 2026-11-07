#!/usr/bin/env python3
"""Measure the Rule 15 code-hygiene criteria over src/.  READ-ONLY -- it changes nothing.

    usage:   python3 tools/auditHygiene.py [--json out.json] [--files]

WHY THIS FILE EXISTS AT ALL.  Five separate hygiene detectors were written against this
codebase and FIVE OF THEM LIED, each by a factor between three and fifty, and each in a way
that looked entirely plausible until a sample was read by eye.  A hygiene count that nobody
can re-derive is worth nothing, so the measurement lives here rather than in a session.
Every calibration below was paid for with a wrong answer; do not simplify one away without
reading a sample of what it changes.

    1. AN ARGUMENT LIST MUST BE SPLIT AT BRACE DEPTH ZERO.  Splitting on every comma reads
       `Array{Tuple{Int64,Float64},1}` as two arguments with no type, and reported 30 % of
       all functions as untyped where the true figure is 1.3 %.
    2. A DOCSTRING HEADER IS COMPARED ON ITS FIRST BACKTICKED NAME, not its last.  Taking
       the last reported 14 header mismatches, every one of them an artifact of a signature
       wrapped over two lines.
    3. A DOCSTRING'S CLOSING TRIPLE QUOTE IS NOT A NEW DEFINITION.  Counting it as one
       reported 4159 definitions with too few blank lines above them; the true figure is 291.
    4. A CONTINUATION LINE IS NOT A BAD INDENT.  A line aligned under an open bracket, or
       under a `using`/`export` list, or a wrapped line of comment prose, is correct where it
       is.  Counting them reported 3282 indent violations where there are, in practice, none.
    5. `Module.` IN A DOCSTRING SIGNATURE IS THE CONVENTION, NOT A DEFECT.  JAC writes
       ``Hfs.Settings(set::Hfs.Settings; ...)`` in the docstring above `function Settings(...)`.
       2017 docstrings differ from their function line for that reason alone.
    6. THE TWO `RacahAlgebra-inc-` FILES ARE OUT OF SCOPE ENTIRELY -- the maintainer's
       decision of 28-Aug-2026.  Their `#` lines are structural to how the sum rules are
       written out.
    7. AN APPENDED `##` IN A DATA TABLE IS A ROW INDEX.  `module-LSjj-inc.jl` alone carries
       20541 of them, numbering the jj-LS coefficient rows against the published tables.
       They are reported separately from ordinary appended comments and are never swept.
"""
import os, re, sys, json, collections

SRC  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")
SKIP = {"module-RacahAlgebra-inc-sumrules.jl", "module-RacahAlgebra-inc-special.jl"}
DATA = {"module-LSjj-inc.jl", "module-LSjj.jl", "module-PeriodicTable.jl",
        "module-Empirical-inc-freefree-gaunt-data.jl", "module-SpinAngularTables-inc-data.jl"}

FUNANY = re.compile(r'^(\s*)function\s+([A-Za-z_][\w.!]*)\s*\(')
FUNCOL = re.compile(r'^function\s+(.+?)\s*$')
GENDEC = re.compile(r'^\s*function\s+[A-Za-z_][\w.!]*\s*end\s*$')
BLOCK  = re.compile(r'^\s*(module|function|for|while|if|let|begin|try|quote|struct|mutable\s+struct|macro)\b')
TRIVIAL = re.compile(r'^(Base\.(show|string|isless|unique|hash|length|getindex|==))')
RETURNED = re.compile(r'is returned|are returned|returns|nothing is returned', re.I)


def docMask(L):
    """-> (inside-docstring flag per line, is-an-opening-delimiter flag per line)."""
    inDoc, doc, opener = False, [False]*len(L), [False]*len(L)
    for i, l in enumerate(L):
        q = l.count('"""')
        if not inDoc and q:
            doc[i] = True
            if q == 1:  inDoc = True;  opener[i] = True
        elif inDoc:
            doc[i] = True
            if q:  inDoc = False
    return doc, opener


def splitArgs(s):
    """Split an argument list at brace, bracket and paren depth zero."""
    out, cur, d = [], "", 0
    for ch in s:
        if   ch in "({[":  d += 1
        elif ch in ")}]":  d -= 1
        if ch in ",;" and d == 0:  out.append(cur);  cur = ""
        else:                      cur += ch
    if cur.strip():  out.append(cur)
    return [a for a in out if a.strip()]


def signature(L, k):
    """The full signature at line k, continued over wrapped lines."""
    sig, t = FUNCOL.match(L[k]).group(1) if FUNCOL.match(L[k]) else L[k].strip(), k
    while sig.count("(") > sig.count(")") and t+1 < len(L):
        t += 1;  sig += " " + L[t].strip()
    return sig


def splitSig(s):
    """-> (bare name, normalised argument list with `Module.` prefixes dropped), or None."""
    p = s.find("(")
    if p < 0:  return None
    name = s[:p].strip()
    name = name.split(".")[-1] if "." in name else name
    args = re.sub(r'\b[A-Za-z_]\w*\.(?=[A-Z])', '', s[p:].split(" where")[0])
    return name, re.sub(r'\s+', '', args)


def audit(path):
    fn = os.path.basename(path)
    L  = open(path, errors="replace").read().split("\n")
    doc, opener = docMask(L)
    r = collections.Counter()
    r["lines"] = len(L)
    sites = collections.defaultdict(list)

    # ---------------------------------------------------------------- comments and whitespace
    depth = 0
    for i, l in enumerate(L):
        s = l.strip()
        if len(l) > 140:
            if doc[i] or s.startswith("#"):  r["longProse"] += 1;  sites["longProse"].append(i+1)
            elif fn in DATA:                 r["longData"]  += 1
            else:                            r["longCode"]  += 1
        if doc[i]:  continue
        if re.match(r'^\s*#\s*$', l):        r["emptyHash"] += 1;  sites["emptyHash"].append(i+1)
        if re.match(r'^\s*#{5,}\s*$', l):    r["hashRule"]  += 1
        elif re.match(r'^\s*#{5,}', l) and re.sub(r'[#\s]', '', l):
            r["textBanner"] += 1;  sites["textBanner"].append(i+1)
        elif re.match(r'^\s*##[^#]', l):     r["standalone2"] += 1;  sites["standalone2"].append(i+1)
        elif re.search(r'\S\s+##[^#]', l) and not s.startswith("#"):
            r["appended2Data" if fn in DATA else "appended2"] += 1
        if re.search(r'==\s*nothing', l):    r["eqNothing"] += 1;  sites["eqNothing"].append(i+1)
        if s and not s.startswith("#") and depth == 0:
            ind = len(l) - len(l.lstrip())
            if ind % 4 and ind < 40:  r["indentOdd"] += 1
        if not s.startswith("#"):
            code = re.sub(r'".*?"', '', l)
            depth = max(0, depth + code.count("(") + code.count("[") - code.count(")") - code.count("]"))

    # ---------------------------------------------------------------- functions
    stack = []
    for i, l in enumerate(L):
        if GENDEC.match(l):  r["genericDecl"] += 1;  continue
        m = FUNANY.match(l)
        if m:
            r["funcs"] += 1
            name = m.group(2)
            r["qualified" if "." in name else "bare"] += 1
            if not m.group(1):                 pass
            elif "function" in stack:          r["nestedFunc"] += 1
            else:  r["indentedTop"] += 1;      sites["indentedTop"].append(i+1)
            args = splitArgs(signature(L, i).split("(", 1)[-1].rsplit(")", 1)[0])
            if args and not all("::" in a or "=" in a or "..." in a for a in args):
                r["untyped"] += 1;  sites["untyped"].append(i+1)
            # docstring immediately above?
            j = i-1
            while j >= 0 and L[j].strip() == "":  j -= 1
            if j >= 0 and L[j].strip().endswith('"""') and doc[j]:
                pass
            elif TRIVIAL.match(name):  r["noDocTrivial"] += 1
            else:                      r["noDoc"] += 1;  sites["noDoc"].append(i+1)
        if doc[i] or l.strip().startswith("#"):  continue
        if BLOCK.match(l) and not re.search(r'\bend\b\s*$', l) and not re.match(r'^\s*\w+\s*\(.*\)\s*=', l):
            stack.append("function" if re.match(r'^\s*function\b', l) else "other")
        elif re.match(r'^\s*end\b', l) and stack:
            stack.pop()

    # ---------------------------------------------------------------- docstring bodies
    i = 0
    while i < len(L):
        if not (opener[i] and L[i].startswith('"""')):  i += 1;  continue
        j = i+1
        while j < len(L) and '"""' not in L[j]:  j += 1
        k = j+1
        while k < len(L) and L[k].strip() == "":  k += 1
        if k < len(L) and FUNCOL.match(L[k]):
            r["funcDocs"] += 1
            body = "\n".join(L[i:j+1])
            if not RETURNED.search(body):  r["noReturnStated"] += 1;  sites["noReturnStated"].append(k+1)
            head = " ".join(L[i:min(i+4, j)])
            cand = next((c for c in re.findall(r'`([^`]+)`', head) if "(" in c), None)
            if cand is None:
                r["noSigInDoc"] += 1;  sites["noSigInDoc"].append(k+1)
            else:
                a, b = splitSig(cand), splitSig(signature(L, k))
                if a and b and a != b:  r["sigMismatch"] += 1;  sites["sigMismatch"].append(k+1)
        i = j+1

    # ---------------------------------------------------------------- blank-line separation
    for i, l in enumerate(L):
        if GENDEC.match(l) or (doc[i] and not opener[i]):  continue
        isDoc = opener[i] and l.startswith('"""')
        isFun = bool(FUNCOL.match(l))
        if not (isDoc or isFun):  continue
        if isFun and i and L[i-1].strip().endswith('"""'):  continue
        n, j = 0, i-1
        while j >= 0 and L[j].strip() == "":  n += 1;  j -= 1
        if j < 0 or re.match(r'^\s*#', L[j]) or L[j].strip().startswith("export"):  continue
        if   n == 0:  r["sep0"] += 1;  sites["sep0"].append(i+1)
        elif n == 1:  r["sep1"] += 1;  sites["sep1"].append(i+1)
        else:         r["sep2plus"] += 1
    return r, sites


ROWS = [("funcs",           "function definitions"),
        ("genericDecl",     "generic `function name end` declarations"),
        ("funcDocs",        "docstrings attached to a function"),
        ("", ""),
        ("standalone2",     "standalone `##` comment lines        -> `#`"),
        ("emptyHash",       "empty `#` lines                      -> a blank line"),
        ("sep1",            "definitions with only 1 blank line   -> 2"),
        ("sep0",            "definitions with no blank line above"),
        ("longProse",       "docstring/comment lines over 140     -> re-flow"),
        ("eqNothing",       "`== nothing` guards                  -> `isnothing`"),
        ("indentedTop",     "top-level functions not at column 1"),
        ("textBanner",      "`#####` separators carrying a sentence"),
        ("", ""),
        ("noDoc",           "functions with no docstring"),
        ("noDocTrivial",    "   of those, trivial Base.* methods (exempt)"),
        ("noReturnStated",  "docstrings not stating the return value"),
        ("sigMismatch",     "docstring signatures disagreeing with the code"),
        ("noSigInDoc",      "docstrings carrying no signature at all"),
        ("untyped",         "functions with an argument lacking a type"),
        ("", ""),
        ("nestedFunc",      "nested helper functions (correctly indented)"),
        ("hashRule",        "pure `#####` separators (kept)"),
        ("longCode",        "code lines over 140 (alignment, not prose)"),
        ("longData",        "data-table lines over 140 (kept)"),
        ("appended2",       "appended `##` on a code line (kept by Rule 15.6)"),
        ("appended2Data",   "   of those, row indices in a data table"),
        ("indentOdd",       "indents not a multiple of 4 (see note 4)"),
        ("qualified",       "`function Module.name(` definitions"),
        ("bare",            "`function name(` definitions")]


def main():
    want  = "--files" in sys.argv
    jpath = sys.argv[sys.argv.index("--json")+1] if "--json" in sys.argv else None
    files = sorted(f for f in os.listdir(SRC) if f.endswith(".jl") and f not in SKIP)
    per, allSites, tot = {}, {}, collections.Counter()
    for f in files:
        r, s = audit(os.path.join(SRC, f))
        per[f] = dict(r);  allSites[f] = {k: v for k, v in s.items()}
        tot.update(r)
    print("src/ hygiene, %d files, %d lines   (the two RacahAlgebra-inc- files are out of scope)"
          % (len(files), tot["lines"]))
    print("-" * 78)
    for key, label in ROWS:
        if not key:  print();  continue
        n = tot[key]
        nf = sum(1 for f in files if per[f].get(key))
        print("   %-52s %6d  in %3d files" % (label, n, nf))
    if want:
        print("\n" + "-" * 78)
        for key, label in ROWS:
            if not key or not tot[key]:  continue
            top = sorted(((per[f].get(key, 0), f) for f in files), reverse=True)[:6]
            print("   %-46s %s" % (label, ", ".join("%s %d" % (f.replace("module-", "").replace(".jl", ""), n)
                                                    for n, f in top if n)))
    if jpath:
        json.dump({"totals": dict(tot), "perFile": per, "sites": allSites}, open(jpath, "w"), indent=0)
        print("\n   per-file counts and line numbers written to %s" % jpath)


main()
