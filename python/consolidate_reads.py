#!/usr/bin/env python3

"""Consolidate per-run FASTA read files into per-sample, per-locus outputs.

Input directory layout:
  <results>/<sample>/<locus>_<run>/reads/
      <sample>_atleast-2.fasta
      <sample>_collapsed_unique.fasta

Output directory layout:
  <output>/<sample>/<locus>/reads/
      <sample>_atleast-2.fasta
      <sample>_collapsed_unique.fasta
"""

from __future__ import annotations

import argparse
import shutil
import sys
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Consolidate <sample>_atleast-2.fasta and "
            "<sample>_collapsed_unique.fasta across runs for each sample/locus."
        )
    )
    parser.add_argument("results_dir", help="Path to source results directory")
    parser.add_argument("output_dir", help="Path to destination output directory")
    return parser.parse_args()


def is_locus_run_dir(name: str) -> tuple[str, str] | None:
    """Return (locus, run) for names matching <locus>_<run>, else None.

    Locus is the text before the first underscore, run is everything after it.
    """
    if "_" not in name:
        return None

    locus, run = name.split("_", 1)
    if not locus or not run:
        return None

    return locus, run


def stream_append(src: Path, dst_handle) -> None:
    """Append file contents from src to an open destination handle."""
    with src.open("rb") as src_handle:
        shutil.copyfileobj(src_handle, dst_handle, length=1024 * 1024)


def warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)


def summarize_sample_loci(sample: str, loci_to_runs: dict[str, list[tuple[str, Path]]]) -> None:
    print(f"Sample: {sample}")
    if not loci_to_runs:
        print("  Loci identified: none")
        return

    print("  Loci identified:")
    for locus in sorted(loci_to_runs):
        runs = sorted(run for run, _ in loci_to_runs[locus])
        print(f"    - {locus}: {len(runs)} run(s) -> {', '.join(runs)}")


def consolidate_sample(
    sample_dir: Path,
    output_root: Path,
) -> tuple[dict[str, dict[str, int]], int]:
    sample = sample_dir.name
    loci_to_runs: dict[str, list[tuple[str, Path]]] = defaultdict(list)

    for child in sorted(sample_dir.iterdir()):
        if not child.is_dir():
            continue

        parsed = is_locus_run_dir(child.name)
        if parsed is None:
            continue

        locus, run = parsed
        loci_to_runs[locus].append((run, child))

    summarize_sample_loci(sample, loci_to_runs)

    sample_summary: dict[str, dict[str, int]] = {}
    missing_files = 0

    for locus in sorted(loci_to_runs):
        output_reads_dir = output_root / sample / locus / "reads"
        output_reads_dir.mkdir(parents=True, exist_ok=True)

        out_atleast = output_reads_dir / f"{sample}_atleast-2.fasta"
        out_collapsed = output_reads_dir / f"{sample}_collapsed_unique.fasta"

        out_atleast.unlink(missing_ok=True)
        out_collapsed.unlink(missing_ok=True)

        appended_atleast = 0
        appended_collapsed = 0
        runs_for_locus = len(loci_to_runs[locus])

        with out_atleast.open("wb") as atleast_handle, out_collapsed.open("wb") as collapsed_handle:
            for run, run_dir in sorted(loci_to_runs[locus], key=lambda item: item[0]):
                reads_dir = run_dir / "reads"

                src_atleast = reads_dir / f"{sample}_atleast-2.fasta"
                src_collapsed = reads_dir / f"{sample}_collapsed_unique.fasta"

                if src_atleast.is_file():
                    stream_append(src_atleast, atleast_handle)
                    appended_atleast += 1
                else:
                    missing_files += 1
                    warn(
                        f"Missing file for sample={sample}, locus={locus}, run={run}: "
                        f"{src_atleast}"
                    )

                if src_collapsed.is_file():
                    stream_append(src_collapsed, collapsed_handle)
                    appended_collapsed += 1
                else:
                    missing_files += 1
                    warn(
                        f"Missing file for sample={sample}, locus={locus}, run={run}: "
                        f"{src_collapsed}"
                    )

        sample_summary[locus] = {
            "runs": runs_for_locus,
            "atleast_files_consolidated": appended_atleast,
            "collapsed_files_consolidated": appended_collapsed,
        }

    print("  Files consolidated:")
    if not sample_summary:
        print("    - none")
    else:
        for locus in sorted(sample_summary):
            stats = sample_summary[locus]
            print(
                "    - "
                f"{locus}: "
                f"{stats['atleast_files_consolidated']}/{stats['runs']} atleast-2, "
                f"{stats['collapsed_files_consolidated']}/{stats['runs']} collapsed_unique"
            )
    print()

    return sample_summary, missing_files


def main() -> int:
    args = parse_args()

    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)

    if not results_dir.is_dir():
        print(f"Error: results_dir does not exist or is not a directory: {results_dir}", file=sys.stderr)
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)

    sample_dirs = sorted(path for path in results_dir.iterdir() if path.is_dir())
    if not sample_dirs:
        print(f"No sample directories found in {results_dir}")
        return 0

    total_samples = 0
    total_loci = 0
    total_missing = 0

    for sample_dir in sample_dirs:
        total_samples += 1
        sample_summary, missing_count = consolidate_sample(sample_dir, output_dir)
        total_loci += len(sample_summary)
        total_missing += missing_count

    print("Consolidation complete.")
    print(f"Samples processed: {total_samples}")
    print(f"Total loci consolidated: {total_loci}")
    print(f"Missing source files: {total_missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
