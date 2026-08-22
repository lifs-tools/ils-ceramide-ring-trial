"""Builders for individual lipidSummarizedQuantities entries.

Pure functions over the dicts returned by ringtrial_sources. Two kinds of entry
are produced: lab-level (one per lab x ceramide x calibration method, carrying
extraction-replicate measurements) and consensus (the published across-lab
values from Suppl_TableS01, which are multi-point).
"""

import ringtrial_terms as terms
from ringtrial_terms import cv

DATASET_FILTER_VALUES = {"All": "All", "Filt": "Without Outliers"}


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
                    metadata, stats, replicates, outlier):
    """One lab x ceramide x calibration-method entry."""
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
        "outlier": cv(*terms.OUTLIER, "true" if outlier else "false"),
    }

    cv_parameters = []
    if metadata["gradientTime"] is not None:
        cv_parameters.append(cv(*terms.GRADIENT_TIME, metadata["gradientTime"]))
    if metadata["sourceTemp"] is not None:
        cv_parameters.append(cv(*terms.SOURCE_TEMPERATURE, metadata["sourceTemp"]))

    measurements = [
        {"quantity": replicate[method],
         "attributes": {"sampleReplicate": cv(*terms.SAMPLE_REPLICATE,
                                              replicate["sample"])}}
        for replicate in replicates if replicate[method] is not None
    ]

    summary = stats[method]
    entry_stats = {"n": len(measurements)}
    if summary["mean"] is not None:
        entry_stats["mean"] = summary["mean"]
    if summary["sd"] is not None:
        entry_stats["sd"] = summary["sd"]

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
        "stats": dict(consensus_row),
    }
