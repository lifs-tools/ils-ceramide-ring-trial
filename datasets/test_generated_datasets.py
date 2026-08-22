"""Invariants asserted against the four generated documents.

These run against in-memory documents, so they do not require a prior
generation run.
"""

import json
import re
import statistics
from pathlib import Path

import jsonschema
import pytest

import generate_datasets as generator
import ringtrial_sources as sources
import ringtrial_terms as terms

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
SCHEMA = json.loads((HERE / "summarized-study-result-schema.json").read_text())

ACCESSION_RE = re.compile(r"^[A-Za-z]+:[A-Za-z0-9]+$")


@pytest.fixture(scope="module")
def documents():
    return generator.build_all(REPO_ROOT)


def test_four_documents_one_per_material(documents):
    assert set(documents) == {"SRM", "hTAG", "DB", "YAA"}


def test_each_document_validates_against_the_schema(documents):
    for document in documents.values():
        jsonschema.validate(document, SCHEMA)


def test_entry_counts(documents):
    for code, document in documents.items():
        quantities = document["lipidSummarizedQuantities"]
        assert len(quantities) == 320, code
        levels = [q["groupingAttributes"]["aggregationLevel"]["value"]
                  for q in quantities]
        assert levels.count("lab") == 312, code
        assert levels.count("consensus") == 8, code


# LabId 20a reported no YAA 4 concentration under either calibration method across
# all four ceramides, so YAA carries eight fewer measurements than the others.
EXPECTED_MEASUREMENTS = {"SRM": 1832, "hTAG": 1832, "DB": 1832, "YAA": 1824}


def test_measurement_count(documents):
    for code, document in documents.items():
        total = sum(len(q.get("individualMeasurements", []))
                    for q in document["lipidSummarizedQuantities"])
        assert total == EXPECTED_MEASUREMENTS[code], code


def test_every_cv_label_is_declared(documents):
    """Regression test for the orphaned-CV bug: an undeclared cvLabel imports
    with a null ControlledVocabulary and never gets a HasCvParent edge."""
    for code, document in documents.items():
        declared = {entry["cvLabel"] for entry in document["controlledVocabularies"]}
        used = set()

        def collect(parameter):
            used.add(parameter["cvLabel"])

        collect(document["studyDesign"])
        for quantity in document["lipidSummarizedQuantities"]:
            collect(quantity["quantityUnit"])
            for parameter in quantity["groupingAttributes"].values():
                collect(parameter)
            for parameter in quantity.get("cvParameters", []):
                collect(parameter)
            for measurement in quantity.get("individualMeasurements", []):
                for parameter in measurement.get("attributes", {}).values():
                    collect(parameter)

        assert used <= declared, f"{code}: undeclared {sorted(used - declared)}"


def test_every_accession_is_well_formed(documents):
    for code, document in documents.items():
        assert ACCESSION_RE.match(document["studyDesign"]["accession"]), \
            (code, document["studyDesign"])
        for quantity in document["lipidSummarizedQuantities"]:
            assert ACCESSION_RE.match(quantity["quantityUnit"]["accession"]), \
                (code, quantity["quantityUnit"])
            for parameter in quantity["groupingAttributes"].values():
                assert ACCESSION_RE.match(parameter["accession"]), (code, parameter)
            for parameter in quantity.get("cvParameters", []):
                assert ACCESSION_RE.match(parameter["accession"]), (code, parameter)
            for measurement in quantity.get("individualMeasurements", []):
                for parameter in measurement.get("attributes", {}).values():
                    assert ACCESSION_RE.match(parameter["accession"]), (code, parameter)


def test_quantity_unit_is_micromole_per_litre_throughout(documents):
    for document in documents.values():
        for quantity in document["lipidSummarizedQuantities"]:
            assert quantity["quantityUnit"]["accession"] == "UO:0010003"


def test_stats_mean_matches_the_measurements(documents):
    """Results table 3 is the mean of Results table 2; the documents must agree."""
    for code, document in documents.items():
        for quantity in document["lipidSummarizedQuantities"]:
            measurements = quantity.get("individualMeasurements")
            if not measurements or quantity["stats"]["n"] < 2:
                continue
            observed = statistics.mean(m["quantity"] for m in measurements)
            assert observed == pytest.approx(quantity["stats"]["mean"], rel=1e-9), code


def test_all_39_labs_appear_in_every_document(documents):
    expected = set(sources.read_lab_metadata(REPO_ROOT / "data" / "definitions"))
    assert len(expected) == 39
    for code, document in documents.items():
        seen = {q["groupingAttributes"]["labId"]["value"]
                for q in document["lipidSummarizedQuantities"]
                if "labId" in q["groupingAttributes"]}
        assert seen == expected, code


def test_lab_id_maps_to_the_expected_lab_num(documents):
    metadata = sources.read_lab_metadata(REPO_ROOT / "data" / "definitions")
    for code, document in documents.items():
        for quantity in document["lipidSummarizedQuantities"]:
            grouping = quantity["groupingAttributes"]
            if "labId" not in grouping:
                continue
            lab_id = grouping["labId"]["value"]
            assert grouping["labNum"]["value"] == metadata[lab_id]["labNum"], code


def test_titles_and_native_ids_are_unique_and_expected(documents):
    titles = {document["title"] for document in documents.values()}
    native_ids = {document["nativeId"] for document in documents.values()}
    assert len(titles) == 4 and len(native_ids) == 4
    assert documents["SRM"]["title"] == terms.MATERIALS["SRM"]["title"]
    assert documents["SRM"]["nativeId"] == "ILS_Ceramide_RingTrial_SRM"


def test_five_labs_carry_vendor_level_instrument_terms(documents):
    vendor_accessions = {accession for accession, _ in terms.VENDOR_MODELS.values()}
    fallback = set()
    for quantity in documents["SRM"]["lipidSummarizedQuantities"]:
        grouping = quantity["groupingAttributes"]
        if "instrumentModel" not in grouping:
            continue
        if grouping["instrumentModel"]["accession"] in vendor_accessions:
            fallback.add(grouping["labId"]["value"])
    assert fallback == {"29a", "29b", "30", "35", "38"}
