#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = []
# ///
"""
Fetch verified BibTeX entries for academic papers via the Crossref API.

Standard library only, so `uv run scripts/fetch_bibtex.py ...` works on a clean
machine with nothing installed.

Two-step workflow, because a title query alone is not trustworthy:

    # 1. SEARCH — look at the candidates and pick the right one
    uv run scripts/fetch_bibtex.py search "Random forests Breiman"

    # 2. FETCH — pull BibTeX for the DOI you actually want
    uv run scripts/fetch_bibtex.py doi 10.1023/A:1010933404324 >> book/references.bib

Why two steps: Crossref indexes software and dataset records alongside papers.
A naive `query.title=Random Forests Breiman&rows=1` returns the *CRAN package*
record (10.32614/cran.package.randomforest, type "dataset"), not Breiman's 2001
paper (10.1023/A:1010933404324). `search` shows you the type, year and authors so
that mistake is visible before it reaches the bibliography.

Options:
    --rows N        candidates to show (default 5)
    --type TYPE     restrict to a Crossref type, e.g. journal-article, book
    --mailto EMAIL  use Crossref's polite pool (faster, more reliable)
    --key KEY       override the generated BibTeX key on `doi`
"""

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

API = "https://api.crossref.org/works"
UA = "MS_LJMR-book-chapter-expander/2.0 (https://github.com/luismadrigal98)"
RATE_LIMIT_SECONDS = 1.0
_last_call = 0.0


def _get(url, accept="application/json", timeout=30, retries=3):
    """GET with a polite user agent, rate limiting and backoff on 429/5xx."""
    global _last_call
    for attempt in range(retries):
        elapsed = time.time() - _last_call
        if elapsed < RATE_LIMIT_SECONDS:
            time.sleep(RATE_LIMIT_SECONDS - elapsed)
        req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": accept})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                _last_call = time.time()
                return resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            _last_call = time.time()
            if e.code in (429, 500, 502, 503, 504) and attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            raise SystemExit(f"error: Crossref returned HTTP {e.code} for {url}")
        except urllib.error.URLError as e:
            _last_call = time.time()
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            raise SystemExit(
                f"error: could not reach Crossref ({e.reason}).\n"
                "       If you are behind a proxy or an allowlist, api.crossref.org\n"
                "       and doi.org both need to be reachable."
            )
    raise SystemExit("error: exhausted retries")


def _authors(item, limit=3):
    names = []
    for a in item.get("author", [])[:limit]:
        fam, given = a.get("family"), a.get("given")
        names.append(f"{fam}, {given}" if fam and given else (fam or given or "?"))
    if len(item.get("author", [])) > limit:
        names.append("et al.")
    return "; ".join(names) if names else "(no author listed)"


def _year(item):
    for field in ("published-print", "published-online", "issued", "created"):
        parts = item.get(field, {}).get("date-parts", [[None]])
        if parts and parts[0] and parts[0][0]:
            return parts[0][0]
    return "n.d."


def cmd_search(args):
    params = {"query.bibliographic": args.query, "rows": str(args.rows)}
    if args.type:
        params["filter"] = f"type:{args.type}"
    if args.mailto:
        params["mailto"] = args.mailto
    data = json.loads(_get(f"{API}?{urllib.parse.urlencode(params)}"))
    items = data.get("message", {}).get("items", [])
    if not items:
        print(f"No results for: {args.query}", file=sys.stderr)
        return 1

    print(f"{len(items)} candidate(s) for {args.query!r} "
          f"— verify before citing:\n", file=sys.stderr)
    for i, it in enumerate(items, 1):
        title = (it.get("title") or ["(untitled)"])[0]
        container = (it.get("container-title") or [""])
        print(f"[{i}] {title}")
        print(f"    {_authors(it)} ({_year(it)})")
        print(f"    type: {it.get('type','?'):<18} in: {container[0] if container else ''}")
        print(f"    DOI:  {it.get('DOI','?')}")
        print()
    print("Pick one, then:  fetch_bibtex.py doi <DOI> >> book/references.bib",
          file=sys.stderr)
    return 0


def _rekey(bibtex, key):
    return re.sub(r"^(\s*@\w+\s*\{)[^,]*,", rf"\g<1>{key},", bibtex, count=1)


def cmd_doi(args):
    doi = args.doi.strip()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if doi.lower().startswith(prefix):
            doi = doi[len(prefix):]
    bib = _get(f"https://doi.org/{urllib.parse.quote(doi)}",
               accept="application/x-bibtex").strip()
    if not bib.startswith("@"):
        raise SystemExit(f"error: no BibTeX returned for DOI {doi}")
    if args.key:
        bib = _rekey(bib, args.key)
    key_match = re.match(r"^\s*@\w+\s*\{([^,]+),", bib)
    if key_match:
        print(f"BibTeX key: {key_match.group(1)}", file=sys.stderr)
    print(bib)
    return 0


def main():
    p = argparse.ArgumentParser(
        description="Fetch verified BibTeX from Crossref (stdlib only).")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("search", help="list candidate works for a query")
    s.add_argument("query", nargs="+")
    s.add_argument("--rows", type=int, default=5)
    s.add_argument("--type", default=None,
                   help="e.g. journal-article, book, book-chapter, proceedings-article")
    s.add_argument("--mailto", default=None)
    s.set_defaults(func=cmd_search)

    d = sub.add_parser("doi", help="emit BibTeX for one DOI")
    d.add_argument("doi")
    d.add_argument("--key", default=None, help="override the BibTeX key")
    d.set_defaults(func=cmd_doi)

    args = p.parse_args()
    if getattr(args, "query", None):
        args.query = " ".join(args.query)
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
