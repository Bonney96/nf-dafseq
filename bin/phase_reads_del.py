#!/usr/bin/env python3

import argparse
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import pysam


VALID_BASES = {"A", "C", "G", "T"}
DAF_BASES = {"C", "G", "R", "Y"}
CONVERTED_BASE_BY_ST = {"CT": "Y", "GA": "R"}
MATCH_CIGAR_OPS = {0, 7, 8}


@dataclass(frozen=True)
class SNPPhaseTarget:
    chrom: str
    pos1: int
    ref_base: str
    alt_base: str
    st_to_use: str


@dataclass(frozen=True)
class DeletionPhaseTarget:
    chrom: str
    start1: int
    end1: int
    overlap_frac: float
    st_to_use: Optional[str] = None

    @property
    def start0(self) -> int:
        return self.start1 - 1

    @property
    def end0(self) -> int:
        return self.end1

    @property
    def interval_length(self) -> int:
        return self.end1 - self.start1 + 1


PhaseTarget = Union[SNPPhaseTarget, DeletionPhaseTarget]


def parse_snp(snp: str) -> Tuple[str, int, str, str]:
    """
    Parse SNP from format: chrom:pos,ref,alt
    Example: chr3:128491529,C,T
    """
    try:
        locus, ref, alt = [part.strip() for part in snp.split(",")]
        chrom, pos_str = [part.strip() for part in locus.split(":")]
        pos = int(pos_str)
    except Exception as exc:
        raise ValueError(
            f"Invalid SNP format: {snp!r}. Expected format: chrom:pos,ref,alt"
        ) from exc

    ref = ref.upper()
    alt = alt.upper()

    if pos < 1:
        raise ValueError(f"SNP position must be 1-based and >= 1, got: {pos}")
    if len(ref) != 1 or ref not in VALID_BASES:
        raise ValueError(f"Invalid reference base: {ref!r}. Use one of A/C/G/T.")
    if len(alt) != 1 or alt not in VALID_BASES:
        raise ValueError(f"Invalid alternate base: {alt!r}. Use one of A/C/G/T.")
    if ref == alt:
        raise ValueError("Reference and alternate base must be different.")

    return chrom, pos, ref, alt


def parse_deletion_interval(deletion: str) -> Tuple[str, int, int]:
    """
    Parse deletion from format: chrom:start-end
    Example: chr3:128489565-128489584
    Coordinates are treated as 1-based inclusive.
    """
    try:
        chrom, coords = [part.strip() for part in deletion.split(":")]
        start_str, end_str = [part.strip() for part in coords.split("-")]
        start1 = int(start_str)
        end1 = int(end_str)
    except Exception as exc:
        raise ValueError(
            f"Invalid deletion format: {deletion!r}. Expected format: chrom:start-end"
        ) from exc

    if start1 < 1 or end1 < 1:
        raise ValueError(
            f"Deletion coordinates must be 1-based and >= 1, got: {start1}-{end1}"
        )
    if end1 < start1:
        raise ValueError(
            f"Deletion end must be >= start, got: {start1}-{end1}"
        )

    interval_length = end1 - start1 + 1
    if interval_length < 3:
        raise ValueError(
            "Deletion interval must span at least 3 nt. "
            f"Received {chrom}:{start1}-{end1} ({interval_length} nt)."
        )

    return chrom, start1, end1


def parse_fraction(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid fraction: {value!r}. Expected a number in the range 0 < value <= 1."
        ) from exc

    if not (0 < parsed <= 1):
        raise argparse.ArgumentTypeError(
            f"Invalid fraction: {value!r}. Expected 0 < value <= 1."
        )

    return parsed


def choose_st_to_use(ref_base: str) -> str:
    return "GA" if ref_base in {"C", "T"} else "CT"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Phase reads by either a SNP or a deletion interval, write a filtered phased "
            "BAM/BAI, and generate non-phased/hap1/hap2 bedGraph and bigWig outputs."
        ),
        epilog=(
            "Examples:\n"
            "  SNP phasing:\n"
            "    phase_reads_v1.py --bam sample.bam --snp chr3:128491529,C,T \\\n"
            "      --results-path results --sample-prefix sample\n"
            "\n"
            "  Deletion phasing:\n"
            "    phase_reads_v1.py --bam sample.bam --del chr3:128489565-128489584 \\\n"
            "      --del-overlap-frac 0.66 --results-path results --sample-prefix sample\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bam", required=True, help="Input BAM path")

    phase_group = parser.add_mutually_exclusive_group(required=True)
    phase_group.add_argument(
        "--snp",
        help='Phase by SNP in format "chrom:pos,ref,alt" (example: chr3:128491529,C,T)',
    )
    phase_group.add_argument(
        "--del",
        dest="deletion",
        help='Phase by deletion interval in format "chrom:start-end" (1-based inclusive; minimum length 3 nt)',
    )

    parser.add_argument(
        "--del-overlap-frac",
        type=parse_fraction,
        default=0.66,
        help=(
            "Minimum fraction of the requested deletion interval that must be covered by "
            "a read's aligned deletion for deletion support. Used only with --del. "
            "Default: 0.66"
        ),
    )
    parser.add_argument(
        "--st-to-use",
        choices=sorted(CONVERTED_BASE_BY_ST),
        default=None,
        help=(
            "Optional strand tag override for phased BAM/bedGraph outputs. "
            "By default, SNP mode preserves the historical SNP-derived choice and "
            "deletion mode auto-selects the predominant phased strand."
        ),
    )
    parser.add_argument("--results-path", required=True, help="Output directory")
    parser.add_argument("--sample-prefix", required=True, help="Prefix for output files")
    parser.add_argument(
        "--chrom-sizes",
        default=None,
        help="Optional chrom sizes file for bedGraphToBigWig. If omitted, generated from BAM header.",
    )
    return parser


def resolve_phase_target(args: argparse.Namespace) -> PhaseTarget:
    if args.snp:
        chrom, pos1, ref_base, alt_base = parse_snp(args.snp)
        st_to_use = args.st_to_use or choose_st_to_use(ref_base)
        return SNPPhaseTarget(
            chrom=chrom,
            pos1=pos1,
            ref_base=ref_base,
            alt_base=alt_base,
            st_to_use=st_to_use,
        )

    chrom, start1, end1 = parse_deletion_interval(args.deletion)
    return DeletionPhaseTarget(
        chrom=chrom,
        start1=start1,
        end1=end1,
        overlap_frac=args.del_overlap_frac,
        st_to_use=args.st_to_use,
    )


def write_chrom_sizes_from_bam_header(bam_path: Path, out_path: Path) -> None:
    with pysam.AlignmentFile(str(bam_path), "rb") as bamfile:
        sq_entries = bamfile.header.to_dict().get("SQ", [])

    with out_path.open("w", encoding="utf-8") as handle:
        for sq in sq_entries:
            handle.write(f"{sq['SN']}\t{sq['LN']}\n")


def write_bedgraph(path: Path, chrom: str, signal_by_pos: Dict[int, float]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for pos0 in sorted(signal_by_pos):
            handle.write(f"{chrom}\t{pos0}\t{pos0 + 1}\t{signal_by_pos[pos0]}\n")


def signal_non_phased(
    base_counts_by_pos: Dict[int, Counter],
) -> Dict[int, float]:
    signal = {}
    for pos0, counts in base_counts_by_pos.items():
        coverage = counts.get("C", 0) + counts.get("G", 0) + counts.get("R", 0) + counts.get("Y", 0)
        if coverage > 0:
            signal[pos0] = (counts.get("R", 0) + counts.get("Y", 0)) / coverage
    return signal


def signal_phased(
    base_counts_by_pos: Dict[int, Counter], st_to_use: str
) -> Dict[int, float]:
    converted_base = CONVERTED_BASE_BY_ST[st_to_use]
    signal = {}
    for pos0, counts in base_counts_by_pos.items():
        coverage = counts.get("C", 0) + counts.get("G", 0) + counts.get(converted_base, 0)
        if coverage > 0:
            signal[pos0] = counts.get(converted_base, 0) / coverage
    return signal


def run_bedgraph_to_bigwig(bedgraph: Path, chrom_sizes: Path, bigwig: Path) -> None:
    cmd = ["bedGraphToBigWig", str(bedgraph), str(chrom_sizes), str(bigwig)]
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "bedGraphToBigWig not found on PATH. In your cluster env run: module load labtools"
        ) from exc


def iter_primary_mapped_reads(bamfile: pysam.AlignmentFile) -> Iterable[pysam.AlignedSegment]:
    # until_eof=True avoids requiring a BAM index on input.
    for read in bamfile.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        yield read


def collect_read_records(
    read: pysam.AlignedSegment, st_value: str, has_du_tag: bool
) -> List[Tuple[str, bool, int, str, str]]:
    seq = read.query_sequence
    if not seq:
        return []

    records = []
    for read_pos, ref_pos, ref_base_at_pos in read.get_aligned_pairs(
        matches_only=True, with_seq=True
    ):
        if read_pos is None or ref_pos is None or ref_base_at_pos is None:
            continue
        records.append(
            (
                st_value,
                has_du_tag,
                ref_pos,
                ref_base_at_pos.upper(),
                seq[read_pos].upper(),
            )
        )
    return records


def assign_haplotype_by_snp(
    per_read_records: List[Tuple[str, bool, int, str, str]],
    target: SNPPhaseTarget,
) -> int:
    target_ref_pos0 = target.pos1 - 1
    for _st_value, _has_du_tag, ref_pos, _ref_base_u, read_base_u in per_read_records:
        if ref_pos != target_ref_pos0:
            continue
        if read_base_u == target.ref_base:
            return 1
        if read_base_u == target.alt_base:
            return 2
    return 0


def deletion_overlap_length(
    read: pysam.AlignedSegment, target: DeletionPhaseTarget
) -> int:
    ref_pos = read.reference_start
    deleted_overlap = 0

    for op, length in read.cigartuples or []:
        if op in MATCH_CIGAR_OPS:
            ref_pos += length
        elif op == 2:
            del_start = ref_pos
            del_end = ref_pos + length
            deleted_overlap += max(
                0, min(del_end, target.end0) - max(del_start, target.start0)
            )
            ref_pos = del_end
        elif op == 3:
            ref_pos += length
        elif op in {1, 4, 5, 6}:
            continue

    return deleted_overlap


def read_spans_deletion_interval(
    read: pysam.AlignedSegment, target: DeletionPhaseTarget
) -> bool:
    if read.reference_start is None or read.reference_end is None:
        return False
    return read.reference_start <= target.start0 and read.reference_end >= target.end0


def assign_haplotype_by_deletion(
    read: pysam.AlignedSegment, target: DeletionPhaseTarget
) -> Tuple[int, int, float, bool]:
    deleted_overlap = deletion_overlap_length(read, target)
    overlap_frac = deleted_overlap / target.interval_length
    spans_interval = read_spans_deletion_interval(read, target)

    if overlap_frac >= target.overlap_frac:
        return 2, deleted_overlap, overlap_frac, spans_interval
    if spans_interval and deleted_overlap == 0:
        return 1, deleted_overlap, overlap_frac, spans_interval
    return 0, deleted_overlap, overlap_frac, spans_interval


def choose_st_to_use_for_deletion(
    bam_path: Path, target: DeletionPhaseTarget
) -> str:
    st_counts = Counter()

    with pysam.AlignmentFile(str(bam_path), "rb") as bamfile:
        for read in iter_primary_mapped_reads(bamfile):
            if read.has_tag("du"):
                continue

            st_value = read.get_tag("st") if read.has_tag("st") else "None"
            if st_value not in CONVERTED_BASE_BY_ST:
                continue

            hp_assigned, _deleted_overlap, _overlap_frac, _spans_interval = (
                assign_haplotype_by_deletion(read, target)
            )
            if hp_assigned in {1, 2}:
                st_counts[st_value] += 1

    if not st_counts:
        raise ValueError(
            "Unable to determine st_to_use automatically for deletion mode. "
            "No phased non-duplicate reads with st tag CT/GA were found. "
            "Provide --st-to-use CT or --st-to-use GA."
        )

    return max(sorted(CONVERTED_BASE_BY_ST), key=lambda st_value: st_counts[st_value])


def update_daf_counts(
    per_read_records: List[Tuple[str, bool, int, str, str]],
    st_to_use: str,
    hp_assigned: int,
    non_phased_counts: Dict[int, Counter],
    hap1_counts: Dict[int, Counter],
    hap2_counts: Dict[int, Counter],
) -> None:
    for rec_st, rec_has_du, rec_ref_pos, rec_ref_base, rec_read_base in per_read_records:
        informative = (rec_st == "CT" and rec_ref_base == "C") or (
            rec_st == "GA" and rec_ref_base == "G"
        )
        if not informative:
            continue
        if rec_read_base not in DAF_BASES:
            continue
        if rec_has_du:
            continue

        non_phased_counts[rec_ref_pos][rec_read_base] += 1

        if rec_st == st_to_use and hp_assigned == 1:
            hap1_counts[rec_ref_pos][rec_read_base] += 1
        elif rec_st == st_to_use and hp_assigned == 2:
            hap2_counts[rec_ref_pos][rec_read_base] += 1


def process_bam(
    bam_path: Path,
    out_bam: Path,
    phase_target: PhaseTarget,
    st_to_use: str,
) -> Dict[str, Any]:
    non_phased_counts = defaultdict(Counter)
    hap1_counts = defaultdict(Counter)
    hap2_counts = defaultdict(Counter)

    phase_counter = {1: 0, 2: 0}
    total_reads = 0
    written_reads = 0
    deletion_ambiguous_reads = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bamfile, pysam.AlignmentFile(
        str(out_bam), "wb", template=bamfile
    ) as out_bam_file:
        for read in iter_primary_mapped_reads(bamfile):
            total_reads += 1

            st_value = read.get_tag("st") if read.has_tag("st") else "None"
            has_du_tag = read.has_tag("du")
            keep_for_bam = (not has_du_tag) and (st_value == st_to_use)
            per_read_records = collect_read_records(read, st_value, has_du_tag)

            if isinstance(phase_target, SNPPhaseTarget):
                hp_assigned = assign_haplotype_by_snp(per_read_records, phase_target)
            else:
                hp_assigned, deleted_overlap, overlap_frac, spans_interval = (
                    assign_haplotype_by_deletion(read, phase_target)
                )
                if (
                    hp_assigned == 0
                    and spans_interval
                    and deleted_overlap > 0
                    and overlap_frac < phase_target.overlap_frac
                ):
                    deletion_ambiguous_reads += 1

            if hp_assigned in {1, 2}:
                read.set_tag("HP", hp_assigned, value_type="i")
                phase_counter[hp_assigned] += 1

            if keep_for_bam:
                out_bam_file.write(read)
                written_reads += 1

            update_daf_counts(
                per_read_records=per_read_records,
                st_to_use=st_to_use,
                hp_assigned=hp_assigned,
                non_phased_counts=non_phased_counts,
                hap1_counts=hap1_counts,
                hap2_counts=hap2_counts,
            )

    pysam.index(str(out_bam))

    return {
        "total_reads": total_reads,
        "written_reads": written_reads,
        "hap1_reads": phase_counter[1],
        "hap2_reads": phase_counter[2],
        "deletion_ambiguous_reads": deletion_ambiguous_reads,
        "non_phased_counts": non_phased_counts,
        "hap1_counts": hap1_counts,
        "hap2_counts": hap2_counts,
    }


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    bam_path = Path(args.bam).resolve()
    results_path = Path(args.results_path).resolve()
    sample_prefix = args.sample_prefix

    if not bam_path.exists():
        raise FileNotFoundError(f"Input BAM does not exist: {bam_path}")

    results_path.mkdir(parents=True, exist_ok=True)
    phase_target = resolve_phase_target(args)

    if args.chrom_sizes:
        chrom_sizes = Path(args.chrom_sizes).resolve()
        if not chrom_sizes.exists():
            raise FileNotFoundError(f"--chrom-sizes file not found: {chrom_sizes}")
    else:
        chrom_sizes = results_path / f"{sample_prefix}.chrom.sizes"
        write_chrom_sizes_from_bam_header(bam_path, chrom_sizes)

    if isinstance(phase_target, DeletionPhaseTarget) and phase_target.st_to_use is None:
        st_to_use = choose_st_to_use_for_deletion(bam_path, phase_target)
        phase_target = DeletionPhaseTarget(
            chrom=phase_target.chrom,
            start1=phase_target.start1,
            end1=phase_target.end1,
            overlap_frac=phase_target.overlap_frac,
            st_to_use=st_to_use,
        )
    else:
        st_to_use = phase_target.st_to_use

    out_bam = results_path / f"{sample_prefix}.filtered_phased.bam"
    out_bai = Path(str(out_bam) + ".bai")
    non_phased_bedgraph = results_path / f"{sample_prefix}.non_phased.bedgraph"
    hap1_bedgraph = results_path / f"{sample_prefix}.hap1.bedgraph"
    hap2_bedgraph = results_path / f"{sample_prefix}.hap2.bedgraph"
    non_phased_bw = results_path / f"{sample_prefix}.non_phased.bw"
    hap1_bw = results_path / f"{sample_prefix}.hap1.bw"
    hap2_bw = results_path / f"{sample_prefix}.hap2.bw"

    results = process_bam(
        bam_path=bam_path,
        out_bam=out_bam,
        phase_target=phase_target,
        st_to_use=st_to_use,
    )

    non_phased_signal = signal_non_phased(results["non_phased_counts"])
    hap1_signal = signal_phased(results["hap1_counts"], st_to_use)
    hap2_signal = signal_phased(results["hap2_counts"], st_to_use)

    write_bedgraph(non_phased_bedgraph, phase_target.chrom, non_phased_signal)
    write_bedgraph(hap1_bedgraph, phase_target.chrom, hap1_signal)
    write_bedgraph(hap2_bedgraph, phase_target.chrom, hap2_signal)

    run_bedgraph_to_bigwig(non_phased_bedgraph, chrom_sizes, non_phased_bw)
    run_bedgraph_to_bigwig(hap1_bedgraph, chrom_sizes, hap1_bw)
    run_bedgraph_to_bigwig(hap2_bedgraph, chrom_sizes, hap2_bw)

    print(f"results_path={results_path}")
    print(f"sample_prefix={sample_prefix}")
    print(f"input_bam={bam_path}")
    if isinstance(phase_target, SNPPhaseTarget):
        print("phase_mode=snp")
        print(
            f"snp={phase_target.chrom}:{phase_target.pos1},{phase_target.ref_base},{phase_target.alt_base}"
        )
        print(f"hap1_reads_with_snp={results['hap1_reads']}")
        print(f"hap2_reads_with_snp={results['hap2_reads']}")
    else:
        print("phase_mode=deletion")
        print(f"deletion={phase_target.chrom}:{phase_target.start1}-{phase_target.end1}")
        print(f"del_overlap_frac={phase_target.overlap_frac}")
        print(f"hap1_reads_without_deletion={results['hap1_reads']}")
        print(f"hap2_reads_with_deletion={results['hap2_reads']}")
        print(f"deletion_ambiguous_reads={results['deletion_ambiguous_reads']}")
    print(f"st_to_use={st_to_use}")
    print(f"total_primary_mapped_reads={results['total_reads']}")
    print(f"reads_written_to_filtered_bam={results['written_reads']}")
    print(f"filtered_phased_bam={out_bam}")
    print(f"filtered_phased_bai={out_bai}")
    print(f"non_phased_bedgraph={non_phased_bedgraph}")
    print(f"hap1_bedgraph={hap1_bedgraph}")
    print(f"hap2_bedgraph={hap2_bedgraph}")
    print(f"non_phased_bigwig={non_phased_bw}")
    print(f"hap1_bigwig={hap1_bw}")
    print(f"hap2_bigwig={hap2_bw}")
    print(f"chrom_sizes={chrom_sizes}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
