#!/usr/bin/env python3
"""Run a shell command; report wall time and the peak RSS of its process tree (RUSAGE_CHILDREN).
usage: l5_runrss.py LOGFILE 'command ...'"""
import resource, subprocess, sys, time
log, cmd = sys.argv[1], sys.argv[2]
t0 = time.time()
rc = subprocess.call(cmd, shell=True, executable='/bin/bash')
ru = resource.getrusage(resource.RUSAGE_CHILDREN)
peak = ru.ru_maxrss / (1e9 if sys.platform == 'darwin' else 1e6)   # bytes on macOS, KiB on Linux
with open(log, 'a') as fh:
    fh.write('%s\trc=%d\twall_s=%.1f\tuser_s=%.1f\tpeak_rss_GB=%.2f\n' % (cmd.replace('\t', ' '), rc, time.time() - t0, ru.ru_utime, peak))
sys.exit(rc)
