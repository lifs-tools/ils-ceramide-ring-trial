"""Unit tests for the ILS ceramide ring trial dataset generator."""

import re

import ringtrial_terms as terms


ACCESSION_RE = re.compile(r"^[A-Za-z]+:[A-Za-z0-9]+$")


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
