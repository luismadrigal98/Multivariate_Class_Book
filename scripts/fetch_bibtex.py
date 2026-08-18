#!/usr/bin/env python3
"""
Fetch BibTeX entries for academic papers via the Crossref API.
Usage: uv run scripts/fetch_bibtex.py "Query title or authors"
"""
import requests
import sys
import urllib.parse

def fetch_bibtex(query):
    url = f"https://api.crossref.org/works?query.title={urllib.parse.quote(query)}&rows=1"
    res = requests.get(url).json()
    items = res.get('message', {}).get('items', [])
    if items:
        doi = items[0]['DOI']
        headers = {'Accept': 'application/x-bibtex'}
        bib = requests.get(f"https://doi.org/{doi}", headers=headers).text
        return bib
    return f"% Could not find {query}"

if __name__ == "__main__":
    if len(sys.argv) > 1:
        query = " ".join(sys.argv[1:])
        print(fetch_bibtex(query))
    else:
        print("Usage: fetch_bibtex.py 'Paper Title'")
