"""Unit tests for the ILS ceramide ring trial dataset generator."""

import re
from pathlib import Path

import pytest

import ringtrial_sources as sources
import ringtrial_terms as terms


ACCESSION_RE = re.compile(r"^[A-Za-z]+:[A-Za-z0-9]+$")

REPO_ROOT = Path(__file__).resolve().parent.parent
XLSX = REPO_ROOT / "Suppl table concentration values.xlsx"
DEFINITIONS = REPO_ROOT / "data" / "definitions"
OUTPUT = REPO_ROOT / "manuscript" / "output"

CERAMIDES = ["Cer 18:1;O2/16:0", "Cer 18:1;O2/18:0",
             "Cer 18:1;O2/24:0", "Cer 18:1;O2/24:1"]


def test_cv_derives_label_from_accession_prefix():
    assert terms.cv("MS:1002582", "QTRAP 6500+") == {
        "accession": "MS:1002582",
        "name": "QTRAP 6500+",
        "cvLabel": "MS",
    }


def test_cv_stringifies_value():
    assert terms.cv("NCIT:C37984", "Laboratory", 1)["value"] == "1"


def test_cv_omits_value_when_none():
    assert "value" not in terms.cv("UO:0010003", "micromole per litre")


def test_declared_vocabularies_cover_every_prefix_used():
    declared = {entry["cvLabel"] for entry in terms.CONTROLLED_VOCABULARIES}
    used = {terms.STUDY_DESIGN["accession"].split(":")[0],
            terms.QUANTITY_UNIT["accession"].split(":")[0]}
    for table in (terms.INSTRUMENT_MODELS, terms.VENDOR_MODELS,
                  terms.ANALYZER_TYPES, terms.CHROMATOGRAPHY, terms.PROTOCOLS):
        used |= {accession.split(":")[0] for accession, _ in table.values()}
    assert used <= declared, f"undeclared vocabularies: {sorted(used - declared)}"


def test_every_accession_is_well_formed():
    accessions = [terms.STUDY_DESIGN["accession"], terms.QUANTITY_UNIT["accession"]]
    for table in (terms.INSTRUMENT_MODELS, terms.VENDOR_MODELS,
                  terms.ANALYZER_TYPES, terms.CHROMATOGRAPHY, terms.PROTOCOLS):
        accessions += [accession for accession, _ in table.values()]
    bad = [a for a in accessions if not ACCESSION_RE.match(a)]
    assert bad == [], f"malformed accessions: {bad}"


def test_quantity_unit_is_micromole_per_litre():
    # UO:0000276 is 'amount per container' and was the bug in the previous generator.
    assert terms.QUANTITY_UNIT["accession"] == "UO:0010003"
    assert terms.QUANTITY_UNIT["name"] == "micromole per litre"


def test_materials_cover_the_four_reference_materials():
    assert set(terms.MATERIALS) == {"SRM", "hTAG", "DB", "YAA"}
    assert terms.MATERIALS["SRM"]["xlsx"] == "NIST SRM"
    assert terms.MATERIALS["SRM"]["display"] == "NIST SRM 1950"
    assert terms.MATERIALS["SRM"]["filename"] == "ils-ceramide-ring-trial-SRM1950.json"


def test_instrument_overrides_name_the_five_disputed_or_missing_labs():
    assert set(terms.INSTRUMENT_OVERRIDES) == {"29a", "29b", "30", "35", "38"}


@pytest.mark.parametrize("raw,expected", [
    ("3", "03"), ("03", "03"), ("10a", "10a"), ("02b", "02b"), ("38", "38"),
])
def test_normalize_lab_id(raw, expected):
    assert sources.normalize_lab_id(raw) == expected


def test_normalize_lab_id_rejects_garbage():
    with pytest.raises(ValueError):
        sources.normalize_lab_id("not-a-lab")


def test_lab_metadata_covers_exactly_the_39_manuscript_labs():
    meta = sources.read_lab_metadata(DEFINITIONS)
    assert len(meta) == 39
    assert meta["02b"]["labNum"] == "1"
    assert meta["02b"]["vendor"] == "Sciex"
    assert meta["02b"]["instrument"] == "QTRAP 6500+"
    assert meta["02b"]["analyzer"] == "QQQ"
    assert meta["02b"]["protocol"] == "SOP"
    assert meta["02b"]["lc"] == "RP"
    assert meta["02b"]["resolution"] == "LowRes"


def test_lab_metadata_strips_trailing_whitespace_in_instrument_names():
    # instruments.csv row for LabId 27 reads "QTRAP 6500 " with a trailing space.
    assert sources.read_lab_metadata(DEFINITIONS)["27"]["instrument"] == "QTRAP 6500"


def test_lab_metadata_maps_na_method_params_to_none():
    meta = sources.read_lab_metadata(DEFINITIONS)
    assert meta["09"]["gradientTime"] is None      # LabProtocolDetails records NA
    assert meta["02b"]["gradientTime"] == "5"
    assert meta["02b"]["sourceTemp"] == "300"


def test_lab_stats_has_one_entry_per_lab_material_ceramide():
    stats = sources.read_lab_stats(XLSX)
    assert len(stats) == 39 * 4 * 4          # 624 groups
    entry = stats[("02b", "SRM", "Cer 18:1;O2/16:0")]
    assert entry["multi-point"]["mean"] == pytest.approx(0.266664952958908)
    assert entry["multi-point"]["sd"] == pytest.approx(0.00682645852132029)


def test_replicates_are_ordered_and_match_the_lab_mean():
    reps = sources.read_replicates(XLSX)[("02b", "SRM", "Cer 18:1;O2/16:0")]
    assert [r["sample"] for r in reps] == [f"SRM {i}" for i in range(1, 7)]
    assert reps[0]["multi-point"] == pytest.approx(0.270695512055809)


def test_every_lab_stat_group_has_replicates():
    stats = sources.read_lab_stats(XLSX)
    reps = sources.read_replicates(XLSX)
    assert set(stats) == set(reps)


def test_outliers_cover_both_calibration_methods():
    flags = sources.read_outliers(OUTPUT)
    assert len(flags) == 624 * 2
    assert sum(1 for k, v in flags.items() if k[3] == "multi-point" and v) == 85
    assert sum(1 for k, v in flags.items() if k[3] == "single-point" and v) == 81


def test_consensus_matches_the_published_table():
    consensus = sources.read_consensus(OUTPUT)
    assert len(consensus) == 4 * 4 * 2       # material x ceramide x {All, Filt}
    row = consensus[("SRM", "Cer 18:1;O2/16:0", "All")]
    assert row["n"] == 39
    assert row["mean"] == pytest.approx(0.2500015504183957)
    filtered = consensus[("SRM", "Cer 18:1;O2/16:0", "Filt")]
    assert filtered["n"] == 35
    assert filtered["mean"] == pytest.approx(0.2438313501584896)
