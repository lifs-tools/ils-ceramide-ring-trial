"""Builders for individual lipidSummarizedQuantities entries.

Pure functions over the dicts returned by ringtrial_sources. Two kinds of entry
are produced: lab-level (one per lab x ceramide x calibration method, carrying
extraction-replicate measurements) and consensus (the published across-lab
values from Suppl_TableS01, which are multi-point).
"""

import ringtrial_terms as terms
from ringtrial_terms import cv

DATASET_FILTER_VALUES = {"All": "All", "Filt": "Without Outliers"}

# Every field the schema allows on a stats object, in reporting order.
STATS_FIELDS = ("n", "mean", "sd", "median", "cv", "rcv",
                "min", "max", "sum", "variance")

# Below this relative difference the workbook and the published CSVs are taken to
# agree. Excel serialises about fifteen significant digits where the CSVs carry
# seventeen, so every matching pair still differs by a few units in the last place.
# Only a real, communicated correction exceeds this -- see correction_factor().
CORRECTION_EPSILON = 1e-6


def stats_block(**values):
    """A stats object carrying every field, null where the source reports nothing.

    Absence and null say different things downstream: an omitted key reads as "this
    study does not use that statistic", an explicit null as "not reported here". Only
    the latter is true, and a uniform shape stops a consumer inferring zero from a
    missing key.
    """
    unknown = sorted(set(values) - set(STATS_FIELDS))
    if unknown:
        raise ValueError(f"unknown stats fields: {unknown}")
    return {field: values.get(field) for field in STATS_FIELDS}


def correction_factor(lab_id, ceramide, method, workbook_mean, published):
    """Scale taking the workbook's concentration to the published one.

    Exactly 1.0 when the two agree within CORRECTION_EPSILON, so the workbook's
    float serialisation never perturbs replicate values by an ulp.

    LabId 26a/26b (LabNum 22) is the one real exception: that laboratory spiked
    three times the intended internal standard amount, so the manuscript multiplies
    its concentrations by three
    (manuscript/manuscript-figures-tables.Rmd:104). The factor is derived rather
    than hard-coded so the generated data cannot drift from the published values.
    """
    if published is None:
        raise ValueError(
            f"LabId {lab_id}: no published concentration for {ceramide} {method}")
    if workbook_mean is None or workbook_mean == 0:
        raise ValueError(
            f"LabId {lab_id}: no workbook mean for {ceramide} {method}; "
            f"cannot scale replicates onto the published concentration")
    factor = published / workbook_mean
    return 1.0 if abs(factor - 1.0) <= CORRECTION_EPSILON else factor


def instrument_term(lab_id, metadata):
    """Model-level PSI-MS term, or the vendor-level term for overridden labs.

    Five laboratories are overridden: two because PSI-MS has no term for their
    instrument, three because instruments.csv and the original lab reports
    disagree. The override is per laboratory, not per model.
    """
    vendor = metadata["vendor"]
    override = terms.INSTRUMENT_OVERRIDES.get(lab_id)
    if override is not None:
        accession, name = terms.VENDOR_MODELS[vendor]
        return cv(accession, name, override)

    key = (vendor, metadata["instrument"])
    if key not in terms.INSTRUMENT_MODELS:
        raise ValueError(
            f"LabId {lab_id}: no PSI-MS term for {key!r}; add it to "
            f"INSTRUMENT_MODELS or list the lab in INSTRUMENT_OVERRIDES")
    accession, name = terms.INSTRUMENT_MODELS[key]
    return cv(accession, name, f"{vendor} {metadata['instrument']}")


def _lookup(table, key, what, lab_id):
    if key not in table:
        raise ValueError(f"LabId {lab_id}: unknown {what}: {key!r}")
    return table[key]


def build_lab_entry(lab_id, material_code, ceramide, method,
                    metadata, stats, replicates, published):
    """One lab x ceramide x calibration-method entry.

    `published` is the manuscript's own concentration and outlier flag for this
    combination. Its concentration is authoritative: the workbook supplies only the
    SD and the per-replicate values, both rescaled onto it by correction_factor().
    """
    material = terms.MATERIALS[material_code]

    analyzer = _lookup(terms.ANALYZER_TYPES, metadata["analyzer"],
                       "mass analyzer type", lab_id)
    chromatography = _lookup(terms.CHROMATOGRAPHY, metadata["lc"],
                             "chromatography mode", lab_id)
    protocol = _lookup(terms.PROTOCOLS, metadata["protocol"], "protocol", lab_id)

    grouping = {
        "aggregationLevel": cv(*terms.AGGREGATION_LEVEL, "lab"),
        "sampleMatrix": cv(*terms.SAMPLE_MATRIX, material["display"]),
        "labNum": cv(*terms.LAB_NUM, metadata["labNum"]),
        "labId": cv(*terms.LAB_ID, lab_id),
        "instrumentModel": instrument_term(lab_id, metadata),
        "massAnalyzerType": cv(*analyzer, metadata["analyzer"]),
        "massAnalyzerResolution": cv(*terms.ANALYZER_RESOLUTION,
                                     metadata["resolution"]),
        "protocol": cv(*protocol, metadata["protocol"]),
        "chromatography": cv(*chromatography, metadata["lc"]),
        "calibrationMethod": cv(*terms.CALIBRATION_METHOD, method),
        "outlier": cv(*terms.OUTLIER, "true" if published["outlier"] else "false"),
    }

    cv_parameters = []
    if metadata["gradientTime"] is not None:
        cv_parameters.append(cv(*terms.GRADIENT_TIME, metadata["gradientTime"]))
    if metadata["sourceTemp"] is not None:
        cv_parameters.append(cv(*terms.SOURCE_TEMPERATURE, metadata["sourceTemp"]))

    summary = stats[method]
    scale = correction_factor(lab_id, ceramide, method,
                              summary["mean"], published["concentration"])

    # Replicates reported as NA are dropped, never coerced to zero: they would
    # otherwise enter both the sum and the denominator of every downstream mean.
    measurements = [
        {"quantity": replicate[method] * scale,
         "attributes": {"sampleReplicate": cv(*terms.SAMPLE_REPLICATE,
                                              replicate["sample"])}}
        for replicate in replicates if replicate[method] is not None
    ]

    entry_stats = stats_block(
        n=len(measurements),
        mean=published["concentration"],
        sd=None if summary["sd"] is None else summary["sd"] * scale,
    )

    entry = {
        "lipids": [ceramide],
        "quantityUnit": dict(terms.QUANTITY_UNIT),
        "groupingAttributes": grouping,
        "stats": entry_stats,
    }
    if cv_parameters:
        entry["cvParameters"] = cv_parameters
    if measurements:
        entry["individualMeasurements"] = measurements
    return entry


def build_consensus_entry(material_code, ceramide, variant, consensus_row):
    """One published across-laboratory consensus entry (multi-point)."""
    material = terms.MATERIALS[material_code]
    return {
        "lipids": [ceramide],
        "quantityUnit": dict(terms.QUANTITY_UNIT),
        "groupingAttributes": {
            "aggregationLevel": cv(*terms.AGGREGATION_LEVEL, "consensus"),
            "sampleMatrix": cv(*terms.SAMPLE_MATRIX, material["display"]),
            "calibrationMethod": cv(*terms.CALIBRATION_METHOD, "multi-point"),
            "datasetFilter": cv(*terms.DATASET_FILTER,
                                DATASET_FILTER_VALUES[variant]),
        },
        "stats": stats_block(**consensus_row),
    }
