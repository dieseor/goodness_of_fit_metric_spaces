#!/usr/bin/env python3
"""Summarize the fast or re-estimated batches used for Table 3."""
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

METHODS = {"fast": "fast_multiplier", "slow": "reestimated"}


def main(method, source, destination):
    expected_method = METHODS[method]
    groups = defaultdict(dict)
    for path in sorted(source.glob("S*_d*_rep*.csv")):
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                key = (row["scenario"], int(row["d"]))
                rep = int(row["rep"])
                # Two re-estimated replications failed and were rerun on their
                # own; only successful runs are kept.
                if row["status"] != "ok":
                    continue
                if row["effective_bootstrap_method"] != expected_method:
                    raise ValueError(f"Unexpected bootstrap method: {path}, {rep}")
                if rep in groups[key]:
                    raise ValueError(f"Duplicate replication {key}, {rep}: {path}")
                groups[key][rep] = row
    if len(groups) != 16 or any(set(rows) != set(range(1, 1001)) for rows in groups.values()):
        raise ValueError("Expected 16 complete groups of 1,000 replications")
    fields = ["scenario", "d", "n", "beta", "ks_rejection_percent",
              "cvm_rejection_percent", "median_elapsed_seconds", "source"]
    if method == "slow":
        fields.insert(4, "replications")
    with destination.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for (scenario, dimension), reps in sorted(groups.items()):
            values = list(reps.values())
            row = dict(scenario=scenario, d=dimension, n=100, beta=0,
                       ks_rejection_percent=100 * sum(r["ks_reject"] == "TRUE" for r in values) / 1000,
                       cvm_rejection_percent=100 * sum(r["cvm_reject"] == "TRUE" for r in values) / 1000,
                       median_elapsed_seconds=statistics.median(float(r["elapsed_seconds"]) for r in values),
                       source=str(source) + "/*.csv")
            if method == "slow":
                row["replications"] = len(values)
            writer.writerow(row)
    print(f"Summarized 16 complete {method} groups (16,000 replications)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("method", choices=sorted(METHODS))
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    main(args.method, args.source, args.destination)
