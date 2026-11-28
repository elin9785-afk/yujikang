"""Compare chromosome pairs with similar GC percentages using the existing summary."""
import argparse
import csv
from itertools import combinations
from math import isfinite
from pathlib import Path


def compare_pairs(rows, max_gc_diff):
    for a, b in combinations(rows, 2):
        low, high = sorted((a, b), key=lambda r: float(r["gc_percent"]))
        gc_diff = float(high["gc_percent"]) - float(low["gc_percent"])
        if gc_diff <= max_gc_diff:
            observed = float(high["mean_gap"]) - float(low["mean_gap"])
            shuffled = float(high["shuffle_mean_gap"]) - float(low["shuffle_mean_gap"])
            yield low["strain"], high["strain"], gc_diff, observed, shuffled


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", nargs="?", type=Path,
                        default=Path("results/strain_summary.csv"))
    parser.add_argument("--max-gc-diff", type=float, default=1.0,
                        help="maximum GC difference in percentage points (default 1)")
    args = parser.parse_args()
    if not isfinite(args.max_gc_diff) or args.max_gc_diff < 0:
        parser.error("GC difference must be finite and nonnegative")
    with args.summary.open(newline="", encoding="utf-8") as source:
        rows = list(csv.DictReader(source))
    print("lower_GC / higher_GC | GC diff (pp) | observed diff | shuffle diff (codons)")
    pairs = sorted(compare_pairs(rows, args.max_gc_diff), key=lambda pair: pair[2])
    for low, high, gc_diff, observed, shuffled in pairs:
        print(f"{low} / {high} | {gc_diff:.3f} | {observed:+.2f} | {shuffled:+.2f}")
    print(f"Matched {len(pairs)} of {len(rows) * (len(rows) - 1) // 2} pairs.")


if __name__ == "__main__":
    main()
