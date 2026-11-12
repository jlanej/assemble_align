#!/usr/bin/env python3
"""Paced GitHub repository search + latest-release lookup (stdlib only).
Usage: github_scout.py queries.tsv repos.txt out.json
queries.tsv: key<TAB>query ; repos.txt: owner/repo per line.
Search API (unauthenticated, 10 req/min) is paced at 7 s; release tags are read from the
github.com/<repo>/releases/latest redirect (no API quota)."""
import sys, json, time, urllib.request, urllib.parse, urllib.error

def api_get(url):
    for attempt in range(4):
        req = urllib.request.Request(url, headers={"Accept": "application/vnd.github+json"})
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.load(r), None
        except urllib.error.HTTPError as e:
            if e.code in (403, 429):
                reset = e.headers.get("X-RateLimit-Reset")
                wait = max(10, int(reset) - int(time.time()) + 2) if reset else 60
                wait = min(wait, 90)
                time.sleep(wait); continue
            return None, f"HTTP {e.code}"
        except Exception as e:
            time.sleep(5); err = str(e)
    return None, "failed after retries"

class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None

def latest_release(full):
    opener = urllib.request.build_opener(NoRedirect)
    req = urllib.request.Request(f"https://github.com/{full}/releases/latest", method="HEAD")
    try:
        r = opener.open(req, timeout=30)
        return {"status": r.status, "tag": None}
    except urllib.error.HTTPError as e:
        loc = e.headers.get("Location", "")
        if e.code in (301, 302) and "/releases/tag/" in loc:
            return {"status": e.code, "tag": urllib.parse.unquote(loc.rsplit("/releases/tag/", 1)[1])}
        return {"status": e.code, "tag": None, "location": loc}
    except Exception as e:
        return {"status": None, "error": str(e)}

def main(qfile, rfile, out):
    res = {"searched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "search": {}, "releases": {}}
    queries = [l.rstrip("\n").split("\t", 1) for l in open(qfile) if l.strip()]
    for key, q in queries:
        url = "https://api.github.com/search/repositories?" + urllib.parse.urlencode({"q": q, "per_page": 5})
        data, err = api_get(url)
        items = []
        if data:
            for it in data.get("items", [])[:5]:
                items.append({k: it.get(k) for k in ("full_name", "html_url", "description", "stargazers_count", "pushed_at", "created_at", "language")})
        res["search"].setdefault(key, []).append({"query": q, "total": (data or {}).get("total_count"), "items": items, "error": err})
        json.dump(res, open(out, "w"), indent=1)
        time.sleep(7)
    repos = [l.strip() for l in open(rfile) if l.strip() and not l.startswith("#")]
    for full in repos:
        res["releases"][full] = latest_release(full)
        time.sleep(1.0)
        json.dump(res, open(out, "w"), indent=1)
    print("done", len(queries), "queries", len(repos), "repos")

if __name__ == "__main__":
    main(*sys.argv[1:4])
