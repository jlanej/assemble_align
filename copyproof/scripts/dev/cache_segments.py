"""Development aid: cache the HMM segmentation of a run (it does not depend on the segment
classifier), so that scripts/dev/try_classifier.py can re-classify in seconds.

  cache_segments.py RUN_DIR
"""
import sys, pickle
from pathlib import Path
import numpy as np
from copyproof.windows import Windows
from copyproof.model import fit_model
from copyproof.cn import CnParams, assessable, emissions, viterbi_star, segments_of_path, state_runs

def main():
    run = Path(sys.argv[1])
    w = Windows.load(str(run / "windows.w1000.npz"))
    model = fit_model(w)
    p = CnParams()
    ok = assessable(w, model, p)
    out = []
    for ci in range(len(w.contigs)):
        sl = w.slice_of(ci)
        ll = emissions(w, model, p, sl, ok)
        path = viterbi_star(ll, p.p_enter, p.p_leave)
        for a, b, d, llr in segments_of_path(path, ll, ok[sl], p):
            if llr < p.min_llr:
                continue
            parts = [(x, y, l) for x, y, l in state_runs(path, a, b, ok[sl], ll) if l >= p.min_llr]
            out.append((ci, a, b, d, llr, parts))
    pickle.dump(out, open(run / "segcache.pkl", "wb"))
    print(run.name, len(out), "segments cached")

if __name__ == "__main__":
    main()
