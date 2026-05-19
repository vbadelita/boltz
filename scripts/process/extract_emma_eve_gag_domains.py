#!/usr/bin/env python3
"""Extract matrix and capsid domain FASTAs from Emma EVE gag CD-search output."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

DEFAULT_ALIGNDATA = Path("emma_eve_gag/aligndata.json")
DEFAULT_INPUT_FASTA = Path("emma_eve_gag/gag_sample_101.fasta")
DEFAULT_OUTPUT_DIR = Path("emma_eve_gag")

OUTPUT_TO_DOMAINS = {
    "matrix_protein": ("Gag_p10", "Gag_MA"),
    "capsid_n_terminal": ("Gag_p24",),
    "capsid_c_terminal": ("Gag_p24_C",),
}

FASTA_LINE_LENGTH = 80


@dataclass(frozen=True)
class DomainHit:
    """A contiguous query-space domain hit extracted from CD-search output."""

    domain_name: str
    start: int
    end: int

    @property
    def length(self) -> int:
        """Return the 1-based inclusive hit length."""

        return self.end - self.start + 1


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for FASTA extraction."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aligndata", type=Path, default=DEFAULT_ALIGNDATA)
    parser.add_argument("--input-fasta", type=Path, default=DEFAULT_INPUT_FASTA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def parse_simple_fasta(path: Path) -> list[tuple[str, str]]:
    """Parse a minimal FASTA file into `(header, sequence)` tuples."""

    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence_parts: list[str] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence_parts)))
            header = line[1:].strip()
            sequence_parts = []
            continue
        sequence_parts.append(line)
    if header is not None:
        records.append((header, "".join(sequence_parts)))
    return records


def load_cdsearch_objects(path: Path) -> list[dict]:
    """Load concatenated JSON objects emitted by the CD-search batch tool."""

    decoder = json.JSONDecoder()
    text = path.read_text()
    objects: list[dict] = []
    index = 0
    while True:
        while index < len(text) and text[index].isspace():
            index += 1
        if index >= len(text):
            return objects
        obj, index = decoder.raw_decode(text, index)
        objects.append(obj)


def extract_title(seq_annot: dict) -> str | None:
    """Return the record title if the CD-search object contains one."""

    for desc in seq_annot.get("desc", []):
        title = desc.get("title")
        if title:
            return title
    return None


def title_to_fasta_header(title: str) -> str:
    """Strip the `Q#n - >` prefix from a CD-search title."""

    _, _, suffix = title.partition(" - >")
    return suffix if suffix else title.removeprefix(">")


def denseg_query_span(denseg: dict) -> tuple[int, int]:
    """Collapse a dense-seg alignment into a contiguous 0-based query span."""

    dim = denseg["dim"]
    segments: list[tuple[int, int]] = []
    for seg_index, seg_len in enumerate(denseg["lens"]):
        query_start = denseg["starts"][seg_index * dim]
        if query_start >= 0:
            segments.append((query_start, query_start + seg_len))
    if not segments:
        raise ValueError("denseg alignment does not contain query residues")
    return segments[0][0], segments[-1][1] - 1


def collect_domain_hits(path: Path) -> dict[str, dict[str, DomainHit]]:
    """Collect the first relevant domain hit per output category for each sequence."""

    domain_to_output = {
        domain_name: output_name
        for output_name, domain_names in OUTPUT_TO_DOMAINS.items()
        for domain_name in domain_names
    }
    hits_by_header: dict[str, dict[str, DomainHit]] = {}
    for obj in load_cdsearch_objects(path):
        seq_annot = obj.get("Seq_annot", {})
        title = extract_title(seq_annot)
        if title is None:
            continue
        fasta_header = title_to_fasta_header(title)
        output_hits: dict[str, DomainHit] = {}
        for align in seq_annot.get("data", {}).get("align", []):
            labels = [
                item["str"]
                for item in align.get("id", [])
                if isinstance(item, dict) and "str" in item
            ]
            matching_domain = next(
                (label for label in labels if label in domain_to_output),
                None,
            )
            if matching_domain is None:
                continue
            output_name = domain_to_output[matching_domain]
            if output_name in output_hits:
                continue
            start, end = denseg_query_span(align["segs"]["denseg"])
            output_hits[output_name] = DomainHit(
                domain_name=matching_domain,
                start=start,
                end=end,
            )
        if output_hits:
            hits_by_header[fasta_header] = output_hits
    return hits_by_header


def wrap_fasta_sequence(sequence: str) -> str:
    """Wrap a sequence to a conventional FASTA line width."""

    return "\n".join(
        sequence[index : index + FASTA_LINE_LENGTH]
        for index in range(0, len(sequence), FASTA_LINE_LENGTH)
    )


def build_output_records(
    fasta_records: list[tuple[str, str]],
    hits_by_header: dict[str, dict[str, DomainHit]],
) -> dict[str, list[tuple[str, str]]]:
    """Slice the source sequences into per-output FASTA records."""

    outputs = {name: [] for name in OUTPUT_TO_DOMAINS}
    for header, sequence in fasta_records:
        for output_name, hit in hits_by_header.get(header, {}).items():
            subsequence = sequence[hit.start : hit.end + 1]
            output_header = (
                f"{header} | {output_name} | {hit.domain_name} | "
                f"residues={hit.start + 1}-{hit.end + 1}"
            )
            outputs[output_name].append((output_header, subsequence))
    return outputs


def write_output_fastas(
    outputs: dict[str, list[tuple[str, str]]],
    output_dir: Path,
) -> dict[str, Path]:
    """Write extracted subsequences to per-domain FASTA files."""

    output_dir.mkdir(parents=True, exist_ok=True)
    written_paths: dict[str, Path] = {}
    for output_name, records in outputs.items():
        output_path = output_dir / f"{output_name}.fasta"
        with output_path.open("w") as handle:
            for header, sequence in records:
                handle.write(f">{header}\n{wrap_fasta_sequence(sequence)}\n")
        written_paths[output_name] = output_path
    return written_paths


def main() -> int:
    """Run the extraction workflow."""

    args = parse_args()
    fasta_records = parse_simple_fasta(args.input_fasta)
    hits_by_header = collect_domain_hits(args.aligndata)
    outputs = build_output_records(fasta_records, hits_by_header)
    written_paths = write_output_fastas(outputs, args.output_dir)
    for output_name, output_path in written_paths.items():
        print(f"{output_name}\t{len(outputs[output_name])}\t{output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
