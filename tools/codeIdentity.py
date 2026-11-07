#!/usr/bin/env python3
"""Prove that a comment- or whitespace-only pass over src/ changed no code.

    usage:   python3 tools/codeIdentity.py --save  <dir>      before the pass
             python3 tools/codeIdentity.py --check <dir>      after it

It writes, for every `src/*.jl`, a canonical form of the file: every comment removed, every
line that is then empty dropped, trailing blanks trimmed, and NOTHING else touched -- leading
whitespace included, so that a change of indent is still visible.  A pass that only rewrites
comments or inserts blank lines reproduces the canonical form byte for byte, which is a
COMPLETE proof for that class of change and costs seconds rather than the eleven minutes a
test-suite run costs.

Two classes of change it is NOT a proof for, and which have to be read by eye instead:
anything that edits code, such as `== nothing` becoming `isnothing(...)`, and any deliberate
re-indentation.

The comment stripper is Julia-aware.  It keeps a `#` that sits inside a string, inside a
character literal, or inside a triple-quoted block, and it understands the `#= ... =#` block
form.  A stripper that simply cuts at the first `#` corrupts every line carrying a format
string such as "#%3i", of which `src/` holds several, and then reports a false difference in
both directions.  The single quote is read as a character literal only where no value can
precede it, since in Julia it is otherwise the adjoint operator.
"""
import os, sys, hashlib

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")


def strip(text):
    out = []
    i, n = 0, len(text)
    inStr = inChar = inTri = False
    blockDepth = 0
    while i < n:
        c = text[i]
        if blockDepth:
            if   text[i:i+2] == "#=":  blockDepth += 1;  i += 2;  continue
            elif text[i:i+2] == "=#":  blockDepth -= 1;  i += 2;  continue
            if c == "\n":  out.append(c)
            i += 1;  continue
        if inTri:
            if text[i:i+3] == '"""':  out.append('"""');  inTri = False;  i += 3;  continue
            out.append(c);  i += 1;  continue
        if inStr:
            if c == "\\" and i+1 < n:  out.append(text[i:i+2]);  i += 2;  continue
            out.append(c)
            if c == '"':  inStr = False
            i += 1;  continue
        if inChar:
            if c == "\\" and i+1 < n:  out.append(text[i:i+2]);  i += 2;  continue
            out.append(c)
            if c == "'":  inChar = False
            i += 1;  continue
        if text[i:i+3] == '"""':  out.append('"""');  inTri = True;  i += 3;  continue
        if text[i:i+2] == "#=":   blockDepth = 1;  i += 2;  continue
        if c == "#":
            while i < n and text[i] != "\n":  i += 1
            continue
        if c == '"':  out.append(c);  inStr = True;  i += 1;  continue
        if c == "'":
            prev = ""
            for t in reversed(out):
                if not t.isspace():  prev = t[-1];  break
            if not (prev.isalnum() or prev in "_)]."):
                out.append(c);  inChar = True;  i += 1;  continue
        out.append(c);  i += 1
    return "\n".join(l.rstrip() for l in "".join(out).split("\n") if l.strip()) + "\n"


def canon():
    return {f: strip(open(os.path.join(SRC, f), errors="replace").read())
            for f in sorted(x for x in os.listdir(SRC) if x.endswith(".jl"))}


def digest(c):
    return (sum(t.count("\n") for t in c.values()),
            hashlib.sha256("".join(c[f] for f in sorted(c)).encode()).hexdigest()[:16])


def main():
    mode, d = sys.argv[1], sys.argv[2]
    cur = canon()
    nl, sha = digest(cur)
    if mode == "--save":
        os.makedirs(d, exist_ok=True)
        for f, t in cur.items():
            open(os.path.join(d, f + ".canon"), "w").write(t)
        print("canonical form of %d files saved:  %d code lines, sha256 %s" % (len(cur), nl, sha))
        return
    bad, missing = [], [f[:-6] for f in os.listdir(d) if f.endswith(".canon")]
    for f, t in cur.items():
        p = os.path.join(d, f + ".canon")
        if not os.path.exists(p):  bad.append((f, "NEW FILE"));  continue
        missing.remove(f)
        old = open(p, errors="replace").read()
        if old == t:  continue
        ol, nw = old.split("\n"), t.split("\n")
        k = next((j for j in range(min(len(ol), len(nw))) if ol[j] != nw[j]), min(len(ol), len(nw)))
        bad.append((f, "%d -> %d code lines, first difference at canonical line %d\n"
                       "        was: %s\n        now: %s"
                    % (len(ol)-1, len(nw)-1, k+1,
                       ol[k][:100] if k < len(ol) else "<end of file>",
                       nw[k][:100] if k < len(nw) else "<end of file>")))
    print("checked %d files" % len(cur))
    if bad or missing:
        for f, why in bad:  print("   DIFFERS  %s: %s" % (f, why))
        for f in missing:   print("   REMOVED  %s" % f)
        sys.exit(1)
    print("   IDENTICAL -- no code changed.  %d code lines, sha256 %s" % (nl, sha))


main()
