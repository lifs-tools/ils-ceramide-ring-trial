# /// script
# requires-python = ">=3.11"
# dependencies = ["openpyxl", "jsonschema"]
# ///
"""Generate LipidCompass SummarizedStudyResult documents from the ILS ceramide
ring trial.

Emits one JSON document per NIST reference material, each carrying per-laboratory
concentrations (single-point and multi-point calibration, with extraction-replicate
measurements) alongside the published across-laboratory consensus values.

Usage, from anywhere:

    uv run datasets/generate_datasets.py
    uv run datasets/generate_datasets.py \
        --out-dir ../lipidcompass-external-data/SummarizedStudyResults

Design: docs/design/ils-ceramide-ring-trial-import/design.md in lifs-tools/lipidcompass.
"""

import argparse
import json
import sys
from pathlib import Path

import jsonschema

import ringtrial_entries as entries
import ringtrial_sources as sources
import ringtrial_terms as terms

HERE = Path(__file__).resolve().parent
DEFAULT_REPO_ROOT = HERE.parent
SCHEMA_PATH = HERE / "summarized-study-result-schema.json"

CERAMIDES = ["Cer 18:1;O2/16:0", "Cer 18:1;O2/18:0",
             "Cer 18:1;O2/24:0", "Cer 18:1;O2/24:1"]
METHODS = ["single-point", "multi-point"]
VARIANTS = ["All", "Filt"]

DESCRIPTION = (
    "Per-laboratory and consensus ceramide concentrations in {display} from the "
    "interlaboratory ceramide ring trial. Reports {labs} laboratory methods "
    "(LabNum 1-34) measuring four ceramides, with single-point and multi-point "
    "calibration, alongside the published consensus values with and without "
    "outliers. Outliers were determined per ceramide and reference material "
    "using Tukey's 1.5 x IQR fences."
)


def build_document(material_code, metadata, stats, replicates, published, consensus):
    """Assemble one SummarizedStudyResult document for a reference material."""
    material = terms.MATERIALS[material_code]
    quantities = []

    for lab_id in sorted(metadata):
        for ceramide in CERAMIDES:
            key = (lab_id, material_code, ceramide)
            if key not in stats:
                raise ValueError(f"missing lab stats for {key}")
            for method in METHODS:
                quantities.append(entries.build_lab_entry(
                    lab_id, material_code, ceramide, method,
                    metadata[lab_id], stats[key], replicates[key],
                    published[key + (method,)],
                ))

    for ceramide in CERAMIDES:
        for variant in VARIANTS:
            quantities.append(entries.build_consensus_entry(
                material_code, ceramide, variant,
                consensus[(material_code, ceramide, variant)],
            ))

    return {
        "nativeId": material["nativeId"],
        "title": material["title"],
        "description": DESCRIPTION.format(display=material["display"],
                                          labs=len(metadata)),
        "controlledVocabularies": terms.CONTROLLED_VOCABULARIES,
        "studyDesign": terms.STUDY_DESIGN,
        "publicationMetadata": terms.PUBLICATION_METADATA,
        "lipidSummarizedQuantities": quantities,
        "visibility": "PUBLIC",
    }


def build_all(repo_root=DEFAULT_REPO_ROOT):
    """Build all four documents, keyed by material code."""
    repo_root = Path(repo_root)
    definitions = repo_root / "data" / "definitions"
    output = repo_root / "manuscript" / "output"
    xlsx = repo_root / "Suppl table concentration values.xlsx"

    metadata = sources.read_lab_metadata(definitions)
    stats = sources.read_lab_stats(xlsx)
    replicates = sources.read_replicates(xlsx)
    published = sources.read_published_concentrations(output)
    consensus = sources.read_consensus(output)

    _report_corrections(stats, published)

    return {code: build_document(code, metadata, stats, replicates,
                                 published, consensus)
            for code in terms.MATERIALS}


def _report_corrections(stats, published):
    """Print any laboratory whose published concentrations rescale the workbook's.

    Silence here would be the dangerous outcome: a correction that appears in the
    manuscript but not in the workbook is exactly the discrepancy this reader exists
    to absorb, so every run says out loud which laboratories it applied one to.
    """
    factors = {}
    for (lab_id, material, ceramide, method), row in published.items():
        summary = stats[(lab_id, material, ceramide)][method]
        factor = entries.correction_factor(lab_id, ceramide, method,
                                           summary["mean"], row["concentration"])
        if factor != 1.0:
            factors.setdefault(lab_id, []).append(factor)
    for lab_id in sorted(factors):
        values = factors[lab_id]
        print(f"LabId {lab_id}: published concentrations are "
              f"{sum(values) / len(values):.4f}x the workbook values "
              f"({len(values)} of 32 combinations); SD and replicates rescaled to match.")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=DEFAULT_REPO_ROOT,
                        help="ring trial repository root (default: parent of this script)")
    parser.add_argument("--out-dir", type=Path, default=HERE,
                        help="where to write the documents (default: alongside this script)")
    args = parser.parse_args(argv)

    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    documents = build_all(args.repo_root)

    # Validate everything before writing anything, so a failure leaves no
    # half-written set of documents behind.
    for code, document in documents.items():
        try:
            jsonschema.validate(document, schema)
        except jsonschema.ValidationError as error:
            print(f"{code}: schema validation failed at "
                  f"{'/'.join(str(p) for p in error.absolute_path)}: {error.message}",
                  file=sys.stderr)
            return 1

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for code, document in documents.items():
        path = args.out_dir / terms.MATERIALS[code]["filename"]
        path.write_text(json.dumps(document, indent=2, ensure_ascii=False) + "\n",
                        encoding="utf-8")
        quantities = document["lipidSummarizedQuantities"]
        levels = [q["groupingAttributes"]["aggregationLevel"]["value"]
                  for q in quantities]
        measurements = sum(len(q.get("individualMeasurements", []))
                           for q in quantities)
        print(f"{path.name}: {len(quantities)} entries "
              f"({levels.count('lab')} lab, {levels.count('consensus')} consensus), "
              f"{measurements} measurements, {path.stat().st_size // 1024} KiB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
