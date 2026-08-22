"""Unit tests for the ILS ceramide ring trial dataset generator."""

import re
from pathlib import Path

import pytest

import ringtrial_entries as entries
import ringtrial_sources as sources
import ringtrial_terms as terms


ACCESSION_RE = re.compile(r"^[A-Za-z]+:[A-Za-z0-9]+$")

REPO_ROOT = Path(__file__).resolve().parent.parent
XLSX = REPO_ROOT / "Suppl table concentration values.xlsx"
DEFINITIONS = REPO_ROOT / "data" / "definitions"
OUTPUT = REPO_ROOT / "manuscript" / "output"


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


def _fixture():
    return (sources.read_lab_metadata(DEFINITIONS),
            sources.read_lab_stats(XLSX),
            sources.read_replicates(XLSX),
            sources.read_outliers(OUTPUT))


def test_instrument_term_uses_the_exact_model_accession():
    metadata = sources.read_lab_metadata(DEFINITIONS)
    term = entries.instrument_term("02b", metadata["02b"])
    assert term["accession"] == "MS:1002582"
    assert term["name"] == "QTRAP 6500+"
    assert term["value"] == "Sciex QTRAP 6500+"


def test_instrument_term_falls_back_to_vendor_for_overridden_labs():
    metadata = sources.read_lab_metadata(DEFINITIONS)
    shimadzu = entries.instrument_term("30", metadata["30"])
    assert shimadzu["accession"] == "MS:1000124"
    assert shimadzu["name"] == "Shimadzu instrument model"
    assert shimadzu["value"] == "Shimadzu LCMS-8060"

    disputed = entries.instrument_term("38", metadata["38"])
    assert disputed["accession"] == "MS:1000126"
    assert disputed["value"] == "Waters Xevo TQ-S"


def test_instrument_term_override_is_per_lab_not_per_model():
    # LabId 38 is overridden, but 24 and 31 share its model and are not.
    metadata = sources.read_lab_metadata(DEFINITIONS)
    assert entries.instrument_term("24", metadata["24"])["accession"] == "MS:1001792"
    assert entries.instrument_term("31", metadata["31"])["accession"] == "MS:1001792"


def test_build_lab_entry_shape():
    metadata, stats, replicates, outliers = _fixture()
    entry = entries.build_lab_entry(
        "02b", "SRM", "Cer 18:1;O2/16:0", "multi-point",
        metadata["02b"],
        stats[("02b", "SRM", "Cer 18:1;O2/16:0")],
        replicates[("02b", "SRM", "Cer 18:1;O2/16:0")],
        outliers[("02b", "SRM", "Cer 18:1;O2/16:0", "multi-point")],
    )
    assert entry["lipids"] == ["Cer 18:1;O2/16:0"]
    assert entry["quantityUnit"]["accession"] == "UO:0010003"

    grouping = entry["groupingAttributes"]
    assert grouping["aggregationLevel"]["value"] == "lab"
    assert grouping["sampleMatrix"]["value"] == "NIST SRM 1950"
    assert grouping["labNum"]["value"] == "1"
    assert grouping["labId"]["value"] == "02b"
    assert grouping["massAnalyzerType"]["accession"] == "MS:1000081"
    assert grouping["massAnalyzerResolution"]["value"] == "LowRes"
    assert grouping["protocol"]["accession"] == "NCIT:C48443"
    assert grouping["chromatography"]["accession"] == "CHMO:0002302"
    assert grouping["calibrationMethod"]["value"] == "multi-point"
    assert grouping["outlier"]["value"] == "false"

    assert entry["stats"]["n"] == 6
    assert entry["stats"]["mean"] == pytest.approx(0.266664952958908)
    assert len(entry["individualMeasurements"]) == 6
    first = entry["individualMeasurements"][0]
    assert first["quantity"] == pytest.approx(0.270695512055809)
    assert first["attributes"]["sampleReplicate"]["value"] == "SRM 1"


def test_build_lab_entry_omits_absent_method_parameters():
    metadata, stats, replicates, outliers = _fixture()
    key = ("09", "SRM", "Cer 18:1;O2/16:0")
    entry = entries.build_lab_entry(
        "09", "SRM", "Cer 18:1;O2/16:0", "multi-point",
        metadata["09"], stats[key], replicates[key],
        outliers[key + ("multi-point",)],
    )
    # LabId 09 records GradientTime NA and runs FIA rather than RP.
    names = {parameter["name"] for parameter in entry["cvParameters"]}
    assert "gradient elution" not in names
    assert entry["groupingAttributes"]["chromatography"]["accession"] == "MS:1000058"


def test_build_consensus_entry_shape():
    consensus = sources.read_consensus(OUTPUT)
    entry = entries.build_consensus_entry(
        "SRM", "Cer 18:1;O2/16:0", "All",
        consensus[("SRM", "Cer 18:1;O2/16:0", "All")])
    grouping = entry["groupingAttributes"]
    assert grouping["aggregationLevel"]["value"] == "consensus"
    assert grouping["datasetFilter"]["value"] == "All"
    assert grouping["calibrationMethod"]["value"] == "multi-point"
    assert "individualMeasurements" not in entry
    assert entry["stats"]["n"] == 39
    assert entry["stats"]["rcv"] == pytest.approx(11.914665756598856)


def test_build_consensus_entry_labels_the_filtered_variant():
    consensus = sources.read_consensus(OUTPUT)
    entry = entries.build_consensus_entry(
        "SRM", "Cer 18:1;O2/16:0", "Filt",
        consensus[("SRM", "Cer 18:1;O2/16:0", "Filt")])
    assert entry["groupingAttributes"]["datasetFilter"]["value"] == "Without Outliers"
