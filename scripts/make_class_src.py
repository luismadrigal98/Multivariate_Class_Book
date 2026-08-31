#!/usr/bin/env python3
"""
Generate class_src/ -- the plain classroom edition -- from R/.

R/ is the source of truth: it keeps the shared helpers (get_data, need,
with_fig, has_pkg) that let the book's figures be rebuilt in batch. The
professor asked the class scripts to use no wrappers at all, so this script
rewrites each one into plain base R:

    need("vegan")                 ->  library(vegan)
    get_data("biodiv")            ->  read.csv("BiodivCountries.csv")
    file.path(data_dir(), "x")    ->  "x"
    with_fig("name", { BODY })    ->  BODY, dedented, with par() saved/restored
    has_pkg("rgl")                ->  requireNamespace("rgl", quietly = TRUE)
    skip_note(what, pkgs)         ->  message("Skipped: ...")
    auc(p, y)                     ->  a plain auc() defined in the script

Run:  python3 scripts/make_class_src.py
"""
import os, re, sys

SRC_DIR, OUT_DIR = "R", "class_src"
SOURCE_LINE = 'if (!exists("get_data")) source(file.path("R", "00_utils.R"))'
MARK = "##  STANDALONE USE (no repository needed)"

# get_data name -> (plain expression, file to ship, extra setup line)
DATA = {
    "iris":           ("iris",   None, 'data(iris)'),
    "mtcars":         ("mtcars", None, 'data(mtcars)'),
    "doubs":          ("doubs",  None, ('ade4', 'data(doubs)          # the Verneaux river-fish tables')),
    "biodiv":         ('read.csv("BiodivCountries.csv", stringsAsFactors = TRUE)', "BiodivCountries.csv", None),
    "biodiv_pc":      ('read.csv("BiodiversityCountriesPCValues.csv", stringsAsFactors = TRUE)', "BiodiversityCountriesPCValues.csv", None),
    "biodiv_gv":      ('read.csv("BiodiversityCountriesBiGv.csv", stringsAsFactors = TRUE)', "BiodiversityCountriesBiGv.csv", None),
    "butterflies":    ('read.csv("ButterfliesQRoo2.csv", stringsAsFactors = TRUE)', "ButterfliesQRoo2.csv", None),
    "taxon":          ('read.csv("taxon.csv", stringsAsFactors = TRUE)', "taxon.csv", None),
    "crawley":        ('read.csv("speciesCrawley3.csv", stringsAsFactors = TRUE)', "speciesCrawley3.csv", None),
    "hanta":          ('read.csv("hanta_virtual.csv")', "hanta_virtual.csv", None),
    "europe":         ('read.csv("CitiesEurope.csv", stringsAsFactors = TRUE)', "CitiesEurope.csv", None),
    "countries_live": ('read.csv("CountriesToLive.csv", stringsAsFactors = TRUE)', "CountriesToLive.csv", None),
    "neotoma":        ('read.csv("NeotomaMorphoEnvir.csv", stringsAsFactors = TRUE)', "NeotomaMorphoEnvir.csv", None),
    "ssa_factor":     ('read.csv("BiodiversityCountriesSSAFactanal2.csv", stringsAsFactors = TRUE)', "BiodiversityCountriesSSAFactanal2.csv", None),
    "insatisf":       ('read.csv("Insatisf2.csv", stringsAsFactors = TRUE)', "Insatisf2.csv", None),
    "limenitis":      ('as.matrix(read.csv("Limenitis_archippus.csv", header = FALSE))', "Limenitis_archippus.csv", None),
    "leukemia":       ('readRDS("leukemiaExpressionSubset.rds")', "leukemiaExpressionSubset.rds", None),
    "pam":            ('read.csv("pam.csv")', "pam.csv", None),
}

AUC_DEF = '''# Area under the ROC curve, from the rank-sum identity (no extra package)
auc <- function(p, y) {
  r <- rank(p)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
'''


def match_paren(text, i):
    """i indexes the '('; return index just past its matching ')'."""
    depth, j, in_str, esc = 0, i, None, False
    while j < len(text):
        c = text[j]
        if in_str:
            if esc:            esc = False
            elif c == "\\":    esc = True
            elif c == in_str:  in_str = None
        elif c in "\"'":       in_str = c
        elif c in "([{":       depth += 1
        elif c in ")]}":
            depth -= 1
            if depth == 0:     return j + 1
        j += 1
    raise ValueError("unbalanced parentheses")


def split_args(s):
    """Split a R argument list on top-level commas."""
    out, depth, cur, in_str, esc = [], 0, "", None, False
    for c in s:
        if in_str:
            cur += c
            if esc:            esc = False
            elif c == "\\":    esc = True
            elif c == in_str:  in_str = None
            continue
        if c in "\"'":         in_str = c; cur += c; continue
        if c in "([{":         depth += 1
        elif c in ")]}":       depth -= 1
        if c == "," and depth == 0:
            out.append(cur); cur = ""
        else:
            cur += c
    if cur.strip(): out.append(cur)
    return [a.strip() for a in out]


def rejoin_else(text):
    """`if (x) foo()\n else bar()` is legal inside { } but not at top level.
    Dedenting a with_fig block can expose exactly that, so put the `else` back
    on the line it belongs to."""
    lines, out = text.split("\n"), []
    for l in lines:
        if re.match(r'^\s*else\b', l) and out:
            k = len(out) - 1
            while k >= 0 and not out[k].strip():
                k -= 1
            if k >= 0:
                out[k] = out[k].rstrip() + " " + l.strip()
                continue
        out.append(l)
    return "\n".join(out)


def dedent(block, indent):
    lines = block.split("\n")
    out = []
    for l in lines:
        if l.startswith(indent):  out.append(l[len(indent):])
        else:                     out.append(l.lstrip() if l.strip() else "")
    return "\n".join(out)


def expand_with_fig(text):
    """Replace every with_fig(...) call with its body, dedented."""
    while True:
        m = re.search(r'^([ \t]*)with_fig\(', text, re.M)
        if not m: return text
        indent, open_i = m.group(1), m.end() - 1
        close = match_paren(text, open_i)
        args = split_args(text[open_i + 1:close - 1])
        # the body is the first argument that is a { } block, else the 2nd positional
        body = None
        for a in args[1:]:
            if a.startswith("{"):
                body = a[1:-1].strip("\n"); break
        if body is None:
            positional = [a for a in args if not re.match(r'^\w+\s*=', a)]
            body = positional[1] if len(positional) > 1 else ""
        body = dedent(body, indent + "  ") if body.startswith(" ") or "\n" in body else body
        body = "\n".join((indent + l if l.strip() else "") for l in body.split("\n"))
        # A wide or tall multi-panel figure needs room. These scripts open no
        # device, so tell the reader how to make one rather than doing it here.
        def num(argname):
            for a in args:
                mm = re.match(r'^%s\s*=\s*([0-9.]+)$' % argname, a)
                if mm: return float(mm.group(1))
            return None
        w, h = num("width"), num("height")
        if (w and w > 10) or (h and h > 8):
            hint = (indent + "# Wide figure: widen the Plot pane, or open a "
                    "sized device first --\n" + indent +
                    "#   dev.new(width = %g, height = %g)" % (w or 7, h or 7))
            body = hint + "\n" + body
        # Layout hygiene. These scripts open no device, so every figure shares
        # one -- and a par(mfrow) left behind by a previous figure (or by a
        # package's own plot method, e.g. raster::plot) would silently split the
        # next one. A block that sets the layout saves and restores it; a block
        # that does not, states the layout it expects.
        if re.search(r'^\s*par\(', body, re.M):
            body = re.sub(r'^(\s*)par\(', r'\1op <- par(', body, count=1, flags=re.M)
            body = body.rstrip() + "\n" + indent + "par(op)"
        else:
            body = indent + "par(mfrow = c(1, 1))\n" + body
        text = text[:m.start()] + body + text[close:]


def convert(path):
    raw = open(path, encoding="utf8").read()
    name = os.path.basename(path)

    # ---- header: keep authorship, drop the standalone block ------------------
    header = raw[:raw.index(MARK)].rstrip()
    header = header[:header.rindex("## ---")].rstrip() if "## ---" in header else header
    body = raw[raw.index(SOURCE_LINE) + len(SOURCE_LINE):]

    # ---- packages -------------------------------------------------------------
    pkgs = []
    def take_need(m):                       # at column 0 -> hoist to the header
        pkgs.extend(re.findall(r'"([^"]+)"', m.group(1)))
        return ""
    body = re.sub(r'^need\((.*?)\)[ \t]*(?:#.*)?$', take_need, body, flags=re.M)

    def inline_need(m):                     # indented -> library() calls in place
        ind = m.group(1)
        return "\n".join(ind + 'library(%s)' % p
                          for p in re.findall(r'"([^"]+)"', m.group(2)))
    body = re.sub(r'^([ \t]+)need\((.*?)\)[ \t]*(?:#.*)?$', inline_need, body, flags=re.M)

    # ---- data -----------------------------------------------------------------
    used, setup, files = [], [], []
    for n in sorted(set(re.findall(r'get_data\("([a-z_0-9]+)"', body))):
        expr, fil, extra = DATA[n]
        used.append(n)
        if fil:   files.append(fil)
        if isinstance(extra, tuple):        # (package to attach, data() line)
            pkgs.append(extra[0]); setup.append(extra[1])
        elif extra:
            setup.append(extra)
        # VAR <- get_data("iris")  ->  drop when VAR is already the dataset name
        body = re.sub(r'^(\s*)%s\s*<-\s*get_data\("%s"\)\s*$' % (re.escape(expr), n),
                      "", body, flags=re.M)
        body = body.replace('get_data("%s")' % n, expr)

    # ---- the remaining helpers -------------------------------------------------
    direct = []
    def strip_data_dir(m):
        direct.append(m.group(1).strip('"'))
        return m.group(1)
    body = re.sub(r'file\.path\(data_dir\(\),\s*("[^"]+")\)', strip_data_dir, body)
    def hp(m):
        ps = re.findall(r'"([^"]+)"', m.group(1))
        return " && ".join('requireNamespace("%s", quietly = TRUE)' % p for p in ps)
    body = re.sub(r'has_pkg\(([^)]*)\)', hp, body)
    def sn(m):
        a = split_args(m.group(1))
        what = a[0].strip()
        pk   = re.findall(r'"([^"]+)"', a[1]) if len(a) > 1 else []
        inst = "  Install with: install.packages(c(%s))" % \
               ", ".join('\\"%s\\"' % p for p in pk) if pk else ""
        if what.startswith('"') and what.endswith('"'):      # a plain literal
            return 'message("Skipped: %s.%s")' % (what[1:-1], inst)
        return 'message("Skipped: ", %s, ".%s")' % (what, inst)
    body = re.sub(r'skip_note\(([^;]*?)\)\s*$', sn, body, flags=re.M)
    body = rejoin_else(expand_with_fig(body))

    needs_auc = re.search(r'(?<![\w.])auc\(', body) is not None

    # files the script OPENS directly by name, rather than through get_data()
    # (iris.data.csv, vars.tif, ...). Only read calls count -- means_iris.csv is
    # written, not read, and must not appear in the "bring these files" list.
    read_calls = r'(?:read\.csv|read\.table|readRDS|rast|file\.exists)\(\s*"([A-Za-z0-9_.-]+)"'
    # a file is optional when the script guards on its existence, either by
    # literal name or through the variable it was assigned to
    guarded = set(re.findall(r'file\.exists\(\s*"([A-Za-z0-9_.-]+)"', body))
    for var in re.findall(r'file\.exists\(\s*([A-Za-z_.][\w.]*)\s*\)', body):
        for m in re.finditer(r'^\s*%s\s*<-\s*"([^"]+)"' % re.escape(var), body, re.M):
            guarded.add(m.group(1))
    for fn in direct + re.findall(read_calls, body):
        if fn not in files and fn + "  (optional)" not in files:
            files.append(fn + "  (optional)" if fn in guarded else fn)

    # ---- assemble --------------------------------------------------------------
    out = [header, "",
           "##  PLAIN CLASSROOM EDITION",
           "##  Generated from R/%s by scripts/make_class_src.py --" % name,
           "##  edit R/%s and regenerate; changes made here will be overwritten." % name,
           "## ============================================================================",
           ""]
    if pkgs:
        out.append("# R packages required")
        for p in dict.fromkeys(pkgs):
            out.append('#install.packages("%s")' % p)
        for p in dict.fromkeys(pkgs):
            out.append("library(%s)" % p)
        out.append("")
    if files:
        out += ["# Working directory -- point this at the folder holding the data files",
                '#   ' + ", ".join(sorted(set(files))),
                'setwd("YOUR/DIRECTORY")', ""]
    if setup:
        out += ["# Built-in data sets used below"] + setup + [""]
    if needs_auc:
        out += [AUC_DEF]
    out.append(body.strip("\n"))
    text = "\n".join(out) + "\n"
    return re.sub(r'\n{4,}', "\n\n\n", text), files


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    manifest = []
    for f in sorted(os.listdir(SRC_DIR)):
        if not f.endswith(".R") or f == "00_utils.R":
            continue
        text, files = convert(os.path.join(SRC_DIR, f))
        open(os.path.join(OUT_DIR, f), "w", encoding="utf8").write(text)
        manifest.append((f, files))
        print("  %-32s %s" % (f, " ".join(files) or "(no data files)"))
    return manifest


if __name__ == "__main__":
    main()
