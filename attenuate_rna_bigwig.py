#!/usr/bin/env python3

from __future__ import annotations

import argparse
import bisect
import gzip
import heapq
import math
import os
import re
import shlex
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

PROJECT_ROOT = Path(__file__).resolve().parent
MPLCONFIGDIR = PROJECT_ROOT / ".mplconfig"
MPLCONFIGDIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyBigWig


ATTRIBUTE_PATTERN = re.compile(r'(\S+)\s+"([^"]*)"')
DEFAULT_ENTRY_BATCH_SIZE = 250_000


@dataclass
class GeneRecord:
    chrom: str
    start: int
    end: int
    strand: str
    gene_id: str
    gene_name: str
    gene_type: str
    max_signal: float = 0.0
    mean_signal: float = 0.0
    mean_nonzero_signal: float = 0.0
    covered_bases: int = 0
    scale_factor: float = 1.0

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def location(self) -> str:
        return f"{self.chrom}:{self.start + 1:,}-{self.end:,}"

    @property
    def description(self) -> str:
        return f"{self.gene_name} ({self.gene_type}, {self.location}, {self.strand})"


class RunLogger:
    def __init__(self, log_path: Path | None):
        self.log_path = log_path
        self.handle = None
        if log_path is not None:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            self.handle = log_path.open("w", buffering=1)

    def log(self, message: str) -> None:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{timestamp}] {message}"
        print(line, flush=True)
        if self.handle is not None:
            self.handle.write(line + "\n")
            self.handle.flush()

    def close(self) -> None:
        if self.handle is not None:
            self.handle.close()
            self.handle = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Attenuate high-expression genes in an RNA bigWig track while keeping the "
            "result as a bigWig. The script also writes gene-level reports, an outlier "
            "threshold plot, and before/after example plots."
        )
    )
    parser.add_argument("--input-bigwig", required=True, help="Input RNA bigWig file.")
    parser.add_argument("--gtf", required=True, help="Gene annotation GTF or GTF.gz file.")
    parser.add_argument("--output-bigwig", required=True, help="Path for the attenuated bigWig.")
    parser.add_argument(
        "--report-dir",
        required=True,
        help="Directory for TSV summaries and PNG plots.",
    )
    parser.add_argument(
        "--threshold",
        required=True,
        type=float,
        help="Target maximum signal. Genes above this max are attenuated.",
    )
    parser.add_argument(
        "--chroms",
        help=(
            "Comma-separated chromosomes to process. Useful for fast subset runs while "
            "tuning thresholds and performance."
        ),
    )
    parser.add_argument(
        "--strand",
        choices=["+", "-", "both", "auto"],
        default="auto",
        help="Gene strand subset to use. 'auto' infers from the bigWig filename when possible.",
    )
    parser.add_argument(
        "--scale-only",
        action="store_true",
        help=(
            "Only attenuate outlier genes. By default, the script also hard-caps any remaining "
            "signal above the threshold outside those genes."
        ),
    )
    parser.add_argument(
        "--example-count",
        type=int,
        default=8,
        help="Number of attenuated genes to render as before/after examples.",
    )
    parser.add_argument(
        "--example-flank",
        type=int,
        default=1_000,
        help="Extra bases to show on each side of example gene plots.",
    )
    parser.add_argument(
        "--top-label-count",
        type=int,
        default=20,
        help="How many top outlier genes to annotate in the summary plot.",
    )
    parser.add_argument(
        "--entry-batch-size",
        type=int,
        default=DEFAULT_ENTRY_BATCH_SIZE,
        help="How many transformed intervals to send to pyBigWig per write batch.",
    )
    parser.add_argument(
        "--log-intervals-every",
        type=int,
        default=DEFAULT_ENTRY_BATCH_SIZE,
        help="Emit in-chromosome write progress after this many source intervals.",
    )
    parser.add_argument(
        "--log-file",
        help="Optional log file. Defaults to REPORT_DIR/run.log.",
    )
    parser.add_argument(
        "--max-zooms",
        type=int,
        default=10,
        help=(
            "Number of zoom levels to build in the output bigWig. Keep this above 0 for "
            "IGV compatibility."
        ),
    )
    return parser.parse_args()


def infer_strand(selection: str, bigwig_path: str) -> set[str]:
    if selection in {"+", "-"}:
        return {selection}
    if selection == "both":
        return {"+", "-"}

    name = Path(bigwig_path).name.lower()
    plus_tokens = ("plus", "_plus", ".plus", "-plus")
    minus_tokens = ("minus", "_minus", ".minus", "-minus")
    if any(token in name for token in plus_tokens):
        return {"+"}
    if any(token in name for token in minus_tokens):
        return {"-"}
    return {"+", "-"}


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def select_chrom_sizes(chrom_sizes: dict[str, int], chrom_arg: str | None) -> list[tuple[str, int]]:
    if not chrom_arg:
        return list(chrom_sizes.items())

    selected: list[tuple[str, int]] = []
    seen: set[str] = set()
    for raw_chrom in chrom_arg.split(","):
        chrom = raw_chrom.strip()
        if not chrom:
            continue
        if chrom not in chrom_sizes:
            raise SystemExit(f"Chromosome '{chrom}' was requested but is not present in the bigWig header.")
        if chrom in seen:
            continue
        selected.append((chrom, chrom_sizes[chrom]))
        seen.add(chrom)

    if not selected:
        raise SystemExit("--chroms was provided but no valid chromosome names were found.")
    return selected


def parse_attributes(raw_attributes: str) -> dict[str, str]:
    return {key: value for key, value in ATTRIBUTE_PATTERN.findall(raw_attributes)}


def parse_genes(
    gtf_path: str,
    chrom_sizes: dict[str, int],
    allowed_strands: set[str],
    allowed_chroms: set[str],
) -> dict[str, list[GeneRecord]]:
    genes_by_chrom: dict[str, list[GeneRecord]] = defaultdict(list)
    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom = fields[0]
            strand = fields[6]
            if chrom not in chrom_sizes or chrom not in allowed_chroms or strand not in allowed_strands:
                continue

            start = max(0, int(fields[3]) - 1)
            end = min(int(fields[4]), chrom_sizes[chrom])
            if end <= start:
                continue

            attrs = parse_attributes(fields[8])
            gene_id = attrs.get("gene_id", f"{chrom}:{start + 1}-{end}")
            gene_name = attrs.get("gene_name", gene_id)
            gene_type = attrs.get("gene_type", attrs.get("gene_biotype", "unknown"))
            genes_by_chrom[chrom].append(
                GeneRecord(
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    gene_id=gene_id,
                    gene_name=gene_name,
                    gene_type=gene_type,
                )
            )

    for chrom_genes in genes_by_chrom.values():
        chrom_genes.sort(key=lambda gene: (gene.start, gene.end, gene.gene_name))
    return genes_by_chrom


def compute_gene_stats(
    bigwig: pyBigWig.pyBigWig,
    genes_by_chrom: dict[str, list[GeneRecord]],
    chrom_order: list[str],
    logger: RunLogger | None = None,
) -> list[GeneRecord]:
    all_genes: list[GeneRecord] = []
    for chrom_index, chrom in enumerate(chrom_order, start=1):
        genes = genes_by_chrom.get(chrom, [])
        if not genes:
            if logger is not None:
                logger.log(f"[stats] chrom {chrom_index}/{len(chrom_order)} {chrom}: no genes in selected annotation subset")
            continue
        intervals = bigwig.intervals(chrom) or []
        if logger is not None:
            logger.log(
                f"[stats] chrom {chrom_index}/{len(chrom_order)} {chrom}: "
                f"genes={len(genes):,}, source_intervals={len(intervals):,}"
            )
        if not intervals:
            all_genes.extend(genes)
            continue

        interval_ends = [end for _, end, _ in intervals]
        for gene in genes:
            interval_index = bisect.bisect_right(interval_ends, gene.start)
            max_signal = 0.0
            weighted_sum = 0.0
            covered_bases = 0

            while interval_index < len(intervals):
                start, end, value = intervals[interval_index]
                if start >= gene.end:
                    break
                overlap_start = max(start, gene.start)
                overlap_end = min(end, gene.end)
                if overlap_end > overlap_start:
                    span = overlap_end - overlap_start
                    numeric_value = float(value)
                    if numeric_value > max_signal:
                        max_signal = numeric_value
                    weighted_sum += numeric_value * span
                    covered_bases += span
                interval_index += 1

            gene.max_signal = max_signal
            gene.covered_bases = covered_bases
            gene.mean_signal = weighted_sum / gene.length if gene.length else 0.0
            gene.mean_nonzero_signal = weighted_sum / covered_bases if covered_bases else 0.0
            all_genes.append(gene)
    return all_genes


def compute_percentiles(values: np.ndarray) -> dict[str, float]:
    if values.size == 0:
        return {}
    percentiles = [50, 90, 95, 99, 99.5, 99.9]
    raw = np.percentile(values, percentiles)
    return {f"p{p}": float(v) for p, v in zip(percentiles, raw, strict=True)}


def assign_scale_factors(genes: Iterable[GeneRecord], threshold: float) -> list[GeneRecord]:
    attenuated = []
    for gene in genes:
        if gene.max_signal > threshold:
            gene.scale_factor = threshold / gene.max_signal
            attenuated.append(gene)
        else:
            gene.scale_factor = 1.0
    attenuated.sort(key=lambda gene: gene.max_signal, reverse=True)
    return attenuated


def active_min_factor(heap: list[float], active_counts: Counter[float]) -> float:
    while heap and active_counts.get(heap[0], 0) == 0:
        heapq.heappop(heap)
    return heap[0] if heap else 1.0


def build_attenuation_segments(
    genes_by_chrom: dict[str, list[GeneRecord]]
) -> dict[str, list[tuple[int, int, float]]]:
    segments_by_chrom: dict[str, list[tuple[int, int, float]]] = {}
    for chrom, genes in genes_by_chrom.items():
        events: dict[int, dict[str, list[float]]] = defaultdict(lambda: {"start": [], "end": []})
        for gene in genes:
            if gene.scale_factor >= 1.0:
                continue
            events[gene.start]["start"].append(gene.scale_factor)
            events[gene.end]["end"].append(gene.scale_factor)

        if not events:
            continue

        active_counts: Counter[float] = Counter()
        heap: list[float] = []
        segments: list[tuple[int, int, float]] = []
        previous_position: int | None = None

        for position in sorted(events):
            factor = active_min_factor(heap, active_counts)
            if previous_position is not None and position > previous_position and factor < 1.0:
                segments.append((previous_position, position, factor))

            for factor in events[position]["end"]:
                active_counts[factor] -= 1
                if active_counts[factor] <= 0:
                    del active_counts[factor]
            for factor in events[position]["start"]:
                active_counts[factor] += 1
                heapq.heappush(heap, factor)

            previous_position = position

        segments_by_chrom[chrom] = segments
    return segments_by_chrom


def iter_transformed_intervals(
    intervals: list[tuple[int, int, float]],
    segments: list[tuple[int, int, float]],
    threshold: float,
    hard_cap: bool,
) -> Iterable[tuple[int, int, float]]:
    if not intervals:
        return

    segment_index = 0
    segment_count = len(segments)

    for start, end, value in intervals:
        if end <= start:
            continue

        numeric_value = float(value)
        unchanged_value = min(numeric_value, threshold) if hard_cap else numeric_value

        while segment_index < segment_count and segments[segment_index][1] <= start:
            segment_index += 1

        if segment_index >= segment_count:
            yield (int(start), int(end), float(unchanged_value))
            continue

        cursor = start
        local_index = segment_index

        while local_index < segment_count and segments[local_index][0] < end:
            seg_start, seg_end, factor = segments[local_index]
            if seg_end <= cursor:
                local_index += 1
                continue

            if seg_start > cursor:
                untouched_end = min(seg_start, end)
                if untouched_end > cursor:
                    yield (int(cursor), int(untouched_end), float(unchanged_value))
                    cursor = untouched_end
                if cursor >= end:
                    break

            overlap_start = max(cursor, seg_start)
            overlap_end = min(end, seg_end)
            if overlap_end > overlap_start:
                adjusted = numeric_value * factor
                if hard_cap:
                    adjusted = min(adjusted, threshold)
                yield (int(overlap_start), int(overlap_end), float(adjusted))
                cursor = overlap_end
                if cursor >= end:
                    break

            if seg_end <= cursor:
                local_index += 1

        if cursor < end:
            yield (int(cursor), int(end), float(unchanged_value))

        segment_index = local_index


class IntervalBatchWriter:
    def __init__(self, output_bw: pyBigWig.pyBigWig, chrom: str, batch_size: int):
        self.output_bw = output_bw
        self.chrom = chrom
        self.batch_size = max(1, batch_size)
        self.starts: list[int] = []
        self.ends: list[int] = []
        self.values: list[float] = []
        self.pending: tuple[int, int, float] | None = None
        self.total_emitted = 0

    def add(self, start: int, end: int, value: float) -> None:
        if end <= start:
            return

        if self.pending is not None:
            pending_start, pending_end, pending_value = self.pending
            if pending_end == start and math.isclose(pending_value, value, rel_tol=1e-9, abs_tol=1e-9):
                self.pending = (pending_start, end, pending_value)
                return
            self._append_pending()

        self.pending = (start, end, value)

    def _append_pending(self) -> None:
        if self.pending is None:
            return
        start, end, value = self.pending
        self.starts.append(start)
        self.ends.append(end)
        self.values.append(value)
        self.pending = None
        if len(self.starts) >= self.batch_size:
            self.flush()

    def flush(self) -> None:
        if self.pending is not None:
            start, end, value = self.pending
            self.starts.append(start)
            self.ends.append(end)
            self.values.append(value)
            self.pending = None

        if not self.starts:
            return

        self.output_bw.addEntries(
            [self.chrom] * len(self.starts),
            self.starts,
            ends=self.ends,
            values=self.values,
            validate=False,
        )
        self.total_emitted += len(self.starts)
        self.starts.clear()
        self.ends.clear()
        self.values.clear()


def write_bigwig(
    input_bigwig_path: str,
    output_bigwig_path: str,
    chrom_sizes: list[tuple[str, int]],
    segments_by_chrom: dict[str, list[tuple[int, int, float]]],
    threshold: float,
    hard_cap: bool,
    entry_batch_size: int,
    log_intervals_every: int,
    max_zooms: int,
    logger: RunLogger | None = None,
) -> tuple[int, int]:
    with pyBigWig.open(input_bigwig_path) as source_bw, pyBigWig.open(output_bigwig_path, "w") as output_bw:
        output_bw.addHeader(chrom_sizes, maxZooms=max_zooms)

        total_source_intervals = 0
        total_emitted_intervals = 0

        for chrom_index, (chrom, _) in enumerate(chrom_sizes, start=1):
            chrom_start_time = time.time()
            intervals = list(source_bw.intervals(chrom) or [])
            if not intervals:
                if logger is not None:
                    logger.log(f"[write] chrom {chrom_index}/{len(chrom_sizes)} {chrom}: no source intervals")
                continue

            total_source_intervals += len(intervals)
            segments = segments_by_chrom.get(chrom, [])
            if logger is not None:
                logger.log(
                    f"[write] chrom {chrom_index}/{len(chrom_sizes)} {chrom}: "
                    f"source_intervals={len(intervals):,}, attenuation_segments={len(segments):,}"
                )

            batch_writer = IntervalBatchWriter(output_bw, chrom, entry_batch_size)
            progress_every = max(1, log_intervals_every)
            segment_index = 0
            segment_count = len(segments)

            for source_interval_index, (start, end, value) in enumerate(intervals, start=1):
                if end <= start:
                    continue

                numeric_value = float(value)
                unchanged_value = min(numeric_value, threshold) if hard_cap else numeric_value

                while segment_index < segment_count and segments[segment_index][1] <= start:
                    segment_index += 1

                if segment_index >= segment_count:
                    batch_writer.add(int(start), int(end), float(unchanged_value))
                else:
                    cursor = start
                    local_index = segment_index

                    while local_index < segment_count and segments[local_index][0] < end:
                        seg_start, seg_end, factor = segments[local_index]
                        if seg_end <= cursor:
                            local_index += 1
                            continue

                        if seg_start > cursor:
                            untouched_end = min(seg_start, end)
                            if untouched_end > cursor:
                                batch_writer.add(int(cursor), int(untouched_end), float(unchanged_value))
                                cursor = untouched_end
                            if cursor >= end:
                                break

                        overlap_start = max(cursor, seg_start)
                        overlap_end = min(end, seg_end)
                        if overlap_end > overlap_start:
                            adjusted = numeric_value * factor
                            if hard_cap:
                                adjusted = min(adjusted, threshold)
                            batch_writer.add(int(overlap_start), int(overlap_end), float(adjusted))
                            cursor = overlap_end
                            if cursor >= end:
                                break

                        if seg_end <= cursor:
                            local_index += 1

                    if cursor < end:
                        batch_writer.add(int(cursor), int(end), float(unchanged_value))

                    segment_index = local_index

                if logger is not None and source_interval_index % progress_every == 0:
                    pending_entries = len(batch_writer.starts) + (1 if batch_writer.pending is not None else 0)
                    logger.log(
                        f"[write] chrom {chrom}: processed_source_intervals="
                        f"{source_interval_index:,}/{len(intervals):,}, "
                        f"pending_output_entries={pending_entries:,}, "
                        f"flushed_output_entries={batch_writer.total_emitted:,}"
                    )

            batch_writer.flush()
            total_emitted_intervals += batch_writer.total_emitted
            if logger is not None:
                logger.log(
                    f"[write] chrom {chrom}: done source_intervals={len(intervals):,}, "
                    f"written_intervals={batch_writer.total_emitted:,}, "
                    f"elapsed={time.time() - chrom_start_time:.1f}s"
                )

    return total_source_intervals, total_emitted_intervals


def write_gene_report(path: Path, genes: list[GeneRecord], threshold: float) -> None:
    with path.open("w") as handle:
        handle.write(
            "\t".join(
                [
                    "gene_id",
                    "gene_name",
                    "gene_type",
                    "chrom",
                    "start_1based",
                    "end_1based",
                    "strand",
                    "gene_length",
                    "covered_bases",
                    "max_signal",
                    "mean_signal",
                    "mean_nonzero_signal",
                    "scale_factor",
                    "attenuated",
                    "threshold",
                ]
            )
            + "\n"
        )
        for gene in sorted(genes, key=lambda item: item.max_signal, reverse=True):
            handle.write(
                "\t".join(
                    [
                        gene.gene_id,
                        gene.gene_name,
                        gene.gene_type,
                        gene.chrom,
                        str(gene.start + 1),
                        str(gene.end),
                        gene.strand,
                        str(gene.length),
                        str(gene.covered_bases),
                        f"{gene.max_signal:.6f}",
                        f"{gene.mean_signal:.6f}",
                        f"{gene.mean_nonzero_signal:.6f}",
                        f"{gene.scale_factor:.8f}",
                        "yes" if gene.scale_factor < 1.0 else "no",
                        f"{threshold:.6f}",
                    ]
                )
                + "\n"
            )


def write_threshold_summary(path: Path, genes: list[GeneRecord], threshold: float) -> None:
    signal_values = np.array([gene.max_signal for gene in genes if gene.max_signal > 0], dtype=float)
    percentiles = compute_percentiles(signal_values)
    with path.open("w") as handle:
        handle.write(f"threshold\t{threshold:.6f}\n")
        handle.write(f"genes_total\t{len(genes)}\n")
        handle.write(f"genes_with_signal\t{signal_values.size}\n")
        if signal_values.size:
            handle.write(f"max_signal\t{signal_values.max():.6f}\n")
            handle.write(f"mean_signal\t{signal_values.mean():.6f}\n")
        for label, value in percentiles.items():
            handle.write(f"{label}\t{value:.6f}\n")


def plot_outlier_summary(
    path: Path,
    genes: list[GeneRecord],
    threshold: float,
    top_label_count: int,
) -> None:
    signal_genes = [gene for gene in genes if gene.max_signal > 0]
    if not signal_genes:
        return

    signal_genes.sort(key=lambda gene: gene.max_signal, reverse=True)
    values_desc = np.array([gene.max_signal for gene in signal_genes], dtype=float)
    values_asc = np.sort(values_desc)
    percentiles = compute_percentiles(values_desc)

    figure, axes = plt.subplots(1, 2, figsize=(16, 7))

    positive_min = max(values_asc[0], 1e-6)
    bins = np.logspace(np.log10(positive_min), np.log10(values_asc.max()), 60)
    axes[0].hist(values_asc, bins=bins, color="#3a7ca5", alpha=0.85, edgecolor="white", linewidth=0.4)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Per-gene max signal")
    axes[0].set_ylabel("Gene count")
    axes[0].set_title("Distribution of gene max signal")
    axes[0].axvline(threshold, color="#d1495b", linestyle="--", linewidth=2, label=f"threshold = {threshold:g}")
    for label in ("p95", "p99", "p99.5", "p99.9"):
        value = percentiles.get(label)
        if value is None:
            continue
        axes[0].axvline(value, linestyle=":", linewidth=1.5, label=f"{label} = {value:.2f}")
    axes[0].legend(fontsize=8)

    ranks = np.arange(1, len(signal_genes) + 1)
    axes[1].scatter(ranks, values_desc, s=12, alpha=0.65, color="#2a9d8f")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Gene rank by max signal")
    axes[1].set_ylabel("Per-gene max signal")
    axes[1].set_title("Top outlier genes")
    axes[1].axhline(threshold, color="#d1495b", linestyle="--", linewidth=2)

    for rank, gene in enumerate(signal_genes[:top_label_count], start=1):
        label = f"{gene.gene_name}\n{gene.gene_type} | {gene.location}"
        axes[1].annotate(
            label,
            xy=(rank, gene.max_signal),
            xytext=(8, 4),
            textcoords="offset points",
            fontsize=8,
            ha="left",
            va="bottom",
        )

    summary_lines = [
        f"genes with signal: {len(signal_genes):,}",
        f"max: {values_desc.max():.2f}",
        f"p99: {percentiles.get('p99', 0.0):.2f}",
        f"p99.5: {percentiles.get('p99.5', 0.0):.2f}",
        f"p99.9: {percentiles.get('p99.9', 0.0):.2f}",
    ]
    axes[1].text(
        0.98,
        0.02,
        "\n".join(summary_lines),
        transform=axes[1].transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "#cccccc"},
    )

    figure.suptitle("RNA bigWig outlier summary", fontsize=16)
    figure.tight_layout()
    figure.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(figure)


def slugify(text: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", text.strip())
    safe = safe.strip("._")
    return safe or "gene"


def stats_to_array(values: list[float | None]) -> np.ndarray:
    return np.array([0.0 if value is None else float(value) for value in values], dtype=float)


def plot_example_genes(
    report_dir: Path,
    original_bigwig_path: str,
    attenuated_bigwig_path: str,
    genes: list[GeneRecord],
    threshold: float,
    flank_bases: int,
    example_count: int,
) -> None:
    if not genes or example_count <= 0:
        return

    examples_dir = report_dir / "attenuated_gene_examples"
    examples_dir.mkdir(parents=True, exist_ok=True)

    with pyBigWig.open(original_bigwig_path) as original_bw, pyBigWig.open(attenuated_bigwig_path) as attenuated_bw:
        chrom_sizes = original_bw.chroms()
        for index, gene in enumerate(genes[:example_count], start=1):
            region_start = max(0, gene.start - flank_bases)
            region_end = min(chrom_sizes[gene.chrom], gene.end + flank_bases)
            region_span = region_end - region_start
            if region_span <= 1:
                continue

            bin_count = int(min(600, max(200, region_span // 5)))
            x_positions = np.linspace(region_start + 1, region_end, num=bin_count, endpoint=False)
            before = stats_to_array(
                original_bw.stats(gene.chrom, region_start, region_end, nBins=bin_count, type="mean", exact=True)
            )
            after = stats_to_array(
                attenuated_bw.stats(gene.chrom, region_start, region_end, nBins=bin_count, type="mean", exact=True)
            )

            original_ymax = max(float(before.max()), threshold, 1e-6) * 1.05
            modified_ymax = max(float(after.max()), threshold, 1e-6) * 1.05

            figure, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
            panel_specs = [
                ("A. Both On Original Scale", True, True, original_ymax),
                ("B. Original Only On Original Scale", True, False, original_ymax),
                ("C. Attenuated Only On Attenuated Scale", False, True, modified_ymax),
            ]

            for axis, (title, show_before, show_after, ymax) in zip(axes, panel_specs, strict=True):
                if show_before:
                    axis.plot(x_positions, before, color="#457b9d", linewidth=1.4, label="original")
                if show_after:
                    axis.plot(x_positions, after, color="#d1495b", linewidth=1.4, label="attenuated")
                axis.axvspan(gene.start + 1, gene.end, color="#f4a261", alpha=0.2, label="gene body")
                axis.axhline(threshold, color="#264653", linestyle="--", linewidth=1.2, label=f"threshold = {threshold:g}")
                axis.set_ylim(0, ymax)
                axis.set_ylabel("Mean signal per bin")
                axis.set_title(title, loc="left", fontsize=11)

            axes[0].legend(loc="upper right", fontsize=8)
            axes[-1].set_xlabel(f"{gene.chrom} position (1-based)")
            axes[0].text(
                0.01,
                0.98,
                (
                    f"{gene.gene_name} | {gene.gene_type} | {gene.location}\n"
                    f"original max {gene.max_signal:.2f} -> scale {gene.scale_factor:.4f}"
                ),
                transform=axes[0].transAxes,
                ha="left",
                va="top",
                fontsize=9,
                bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "#dddddd"},
            )
            figure.suptitle(gene.description, fontsize=13)
            figure.tight_layout(rect=(0, 0, 1, 0.97))
            output_name = f"{index:02d}_{slugify(gene.gene_name)}.png"
            figure.savefig(examples_dir / output_name, dpi=200, bbox_inches="tight")
            plt.close(figure)


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def main() -> int:
    args = parse_args()
    if args.threshold <= 0:
        raise SystemExit("--threshold must be greater than 0.")
    if args.entry_batch_size <= 0:
        raise SystemExit("--entry-batch-size must be greater than 0.")
    if args.log_intervals_every <= 0:
        raise SystemExit("--log-intervals-every must be greater than 0.")
    if args.max_zooms < 0:
        raise SystemExit("--max-zooms must be 0 or greater.")

    report_dir = Path(args.report_dir)
    report_dir.mkdir(parents=True, exist_ok=True)
    output_bigwig_path = Path(args.output_bigwig)
    ensure_parent(output_bigwig_path)
    log_path = Path(args.log_file) if args.log_file else report_dir / "run.log"

    allowed_strands = infer_strand(args.strand, args.input_bigwig)
    hard_cap = not args.scale_only
    logger = RunLogger(log_path)

    try:
        logger.log(f"[start] command={shlex.join(sys.argv)}")
        logger.log(
            f"[start] input_bigwig={args.input_bigwig} threshold={args.threshold:g} "
            f"strand={''.join(sorted(allowed_strands))} hard_cap={'yes' if hard_cap else 'no'} "
            f"max_zooms={args.max_zooms}"
        )

        with pyBigWig.open(args.input_bigwig) as input_bw:
            all_chrom_sizes = input_bw.chroms()
            selected_chrom_sizes = select_chrom_sizes(all_chrom_sizes, args.chroms)
            selected_chroms = [chrom for chrom, _ in selected_chrom_sizes]
            logger.log(
                f"[start] selected_chromosomes={','.join(selected_chroms)} "
                f"count={len(selected_chroms):,}"
            )
            genes_by_chrom = parse_genes(
                args.gtf,
                all_chrom_sizes,
                allowed_strands,
                set(selected_chroms),
            )
            genes = compute_gene_stats(
                input_bw,
                genes_by_chrom,
                chrom_order=selected_chroms,
                logger=logger,
            )

        if not genes:
            raise SystemExit("No genes from the GTF matched the requested chromosomes/strand.")

        attenuated_genes = assign_scale_factors(genes, args.threshold)
        segments_by_chrom = build_attenuation_segments(genes_by_chrom)
        logger.log(
            f"[attenuation] genes_analyzed={len(genes):,} genes_attenuated={len(attenuated_genes):,}"
        )

        logger.log("[reports] writing gene-level tables and threshold plot")
        write_gene_report(report_dir / "gene_signal_report.tsv", genes, args.threshold)
        write_threshold_summary(report_dir / "threshold_summary.tsv", genes, args.threshold)
        plot_outlier_summary(
            report_dir / "outlier_threshold_plot.png",
            genes=genes,
            threshold=args.threshold,
            top_label_count=args.top_label_count,
        )

        logger.log("[write] starting bigWig emission from original intervals only")
        total_source_intervals, total_emitted_intervals = write_bigwig(
            input_bigwig_path=args.input_bigwig,
            output_bigwig_path=str(output_bigwig_path),
            chrom_sizes=selected_chrom_sizes,
            segments_by_chrom=segments_by_chrom,
            threshold=args.threshold,
            hard_cap=hard_cap,
            entry_batch_size=args.entry_batch_size,
            log_intervals_every=args.log_intervals_every,
            max_zooms=args.max_zooms,
            logger=logger,
        )
        logger.log(
            f"[write] completed source_intervals={total_source_intervals:,} "
            f"written_intervals={total_emitted_intervals:,}"
        )

        logger.log("[plots] writing before/after examples for attenuated genes")
        plot_example_genes(
            report_dir=report_dir,
            original_bigwig_path=args.input_bigwig,
            attenuated_bigwig_path=str(output_bigwig_path),
            genes=attenuated_genes,
            threshold=args.threshold,
            flank_bases=args.example_flank,
            example_count=args.example_count,
        )

        signal_values = np.array([gene.max_signal for gene in genes if gene.max_signal > 0], dtype=float)
        percentiles = compute_percentiles(signal_values)
        logger.log(f"[done] output_bigwig={output_bigwig_path}")
        logger.log(f"[done] report_dir={report_dir}")
        logger.log(f"[done] genes_analyzed={len(genes):,} genes_attenuated={len(attenuated_genes):,}")
        if signal_values.size:
            logger.log(f"[done] observed_max_signal={signal_values.max():.4f}")
        for label in ("p95", "p99", "p99.5", "p99.9"):
            if label in percentiles:
                logger.log(f"[done] {label}={percentiles[label]:.4f}")
        logger.log(f"[done] hard_cap_enabled={'yes' if hard_cap else 'no'}")
        logger.log(f"[done] log_file={log_path}")
        return 0
    finally:
        logger.close()


if __name__ == "__main__":
    sys.exit(main())
