#!/usr/bin/env python3
"""Summarize the matched fast batches used for Table 3."""
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def main(source, destination):
    groups = defaultdict(dict)
    for path in sorted(source.glob("S*_d*_rep*.csv")):
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                key = (row["scenario"], int(row["d"]))
                rep = int(row["rep"])
                if rep in groups[key]:
                    raise ValueError(f"Duplicate replication {key}, {rep}: {path}")
                if row["status"] != "ok" or row["effective_bootstrap_method"] != "fast_multiplier":
                    raise ValueError(f"Unsuccessful or non-fast replication: {path}, {rep}")
                groups[key][rep] = row
    if len(groups) != 16 or any(set(rows) != set(range(1, 1001)) for rows in groups.values()):
        raise ValueError("Expected 16 complete groups of 1,000 replications")
    fields = ("scenario", "d", "n", "beta", "ks_rejection_percent",
              "cvm_rejection_percent", "median_elapsed_seconds", "source")
    with destination.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for (scenario, dimension), reps in sorted(groups.items()):
            values = list(reps.values())
            writer.writerow(dict(scenario=scenario, d=dimension, n=100, beta=0,
                                 ks_rejection_percent=100 * sum(r["ks_reject"] == "TRUE" for r in values) / 1000,
                                 cvm_rejection_percent=100 * sum(r["cvm_reject"] == "TRUE" for r in values) / 1000,
                                 median_elapsed_seconds=statistics.median(float(r["elapsed_seconds"]) for r in values),
                                 source=str(source) + "/*.csv"))
    print("Summarized 16 complete fast groups (16,000 replications)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    main(args.source, args.destination)
