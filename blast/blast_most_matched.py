#!/usr/bin/env python
import sys

def increment_counts(qseqid, sseqid, qstart, qend, bonus):
    global current_query, current_counts
    assert qseqid == current_query
    if sseqid not in current_counts:
        current_counts[sseqid] = set(), bonus
    else:
        counts, old_bonus = current_counts[sseqid]
        new_bonus = []
        for old, new in zip(old_bonus, bonus):
            if old==new:
                new_bonus.append(old)
            else:
                new_bonus.append(None)
        current_counts[sseqid] =  counts, new_bonus
    current_counts[sseqid][0].update(range(qstart - 1, qend))

def report_old():
    #print current_query, current_counts
    top, bases, bonus = sorted((len(c), k, b) for k, (c, b) in current_counts.items())[0]
    print current_query, top, bases, bonus

current_query = None
current_counts = dict()
bonus_cols = None
for line in sys.stdin:
    parts = line.rstrip("\n").split("\t")
    assert len(parts) > 12
    qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = parts[:12]
    if qseqid != current_query:
        if current_query:
            report_old()
        current_query = qseqid
        current_counts = dict()
        bonus_cols = None
    increment_counts(qseqid, sseqid, int(qstart), int(qend), parts[12:])
if current_query:
    report_old()
