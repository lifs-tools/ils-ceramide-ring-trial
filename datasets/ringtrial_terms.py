"""Controlled-vocabulary term tables for the ILS ceramide ring trial datasets.

Every accession here was resolved against OLS4 (https://www.ebi.ac.uk/ols4) and its
label verified to match. Provenance and rationale: section 6 of
docs/design/ils-ceramide-ring-trial-import/design.md in lifs-tools/lipidcompass.

Do not substitute accessions. If a term looks wrong, raise it rather than changing it.
"""

CONTROLLED_VOCABULARIES = [
    {"name": "Ontology for Biomedical Investigations",
     "uri": "http://purl.obolibrary.org/obo/obi.owl",
     "version": "2026-07-27", "cvLabel": "OBI"},
    {"name": "PSI Mass Spectrometry Ontology",
     "uri": "http://purl.obolibrary.org/obo/ms.owl",
     "version": "4.1.257", "cvLabel": "MS"},
    {"name": "Chemical Methods Ontology",
     "uri": "http://purl.obolibrary.org/obo/chmo.owl",
     "version": "2026-05-28", "cvLabel": "CHMO"},
    {"name": "National Cancer Institute Thesaurus",
     "uri": "http://purl.obolibrary.org/obo/ncit.owl",
     "version": "26.02d", "cvLabel": "NCIT"},
    {"name": "Statistics Ontology",
     "uri": "http://purl.obolibrary.org/obo/stato.owl",
     "version": "2026-04-20", "cvLabel": "STATO"},
    {"name": "Units of Measurement Ontology",
     "uri": "http://purl.obolibrary.org/obo/uo.owl",
     "version": "2026-07-31", "cvLabel": "UO"},
]


def cv(accession, name, value=None):
    """Build a CvParameter dict. cvLabel is always the accession prefix."""
    term = {"accession": accession, "name": name}
    if value is not None:
        term["value"] = str(value)
    term["cvLabel"] = accession.split(":", 1)[0]
    return term


STUDY_DESIGN = cv("OBI:0003658", "ring trial design", "ring trial design")
QUANTITY_UNIT = cv("UO:0010003", "micromole per litre")

# Fixed grouping-attribute terms. Values are supplied per entry.
AGGREGATION_LEVEL = ("NCIT:C95130", "Interlaboratory")
SAMPLE_MATRIX = ("CHMO:0002879", "certified reference material")
LAB_NUM = ("NCIT:C37984", "Laboratory")
LAB_ID = ("NCIT:C93857", "Laboratory Identifier")
ANALYZER_RESOLUTION = ("MS:1000011", "mass resolution")
CALIBRATION_METHOD = ("NCIT:C69187", "Calibration")
OUTLIER = ("STATO:0000036", "outlier")
DATASET_FILTER = ("OBI:0001755", "selection criterion")
SAMPLE_REPLICATE = ("OBI:0000747", "material sample")
GRADIENT_TIME = ("CHMO:0002877", "gradient elution")
SOURCE_TEMPERATURE = ("MS:1002041", "source temperature")

MATERIALS = {
    "SRM": {
        "xlsx": "NIST SRM",
        "display": "NIST SRM 1950",
        "title": "Ceramide Ring Trial — NIST SRM 1950",
        "nativeId": "ILS_Ceramide_RingTrial_SRM",
        "filename": "ils-ceramide-ring-trial-SRM1950.json",
    },
    "hTAG": {
        "xlsx": "NIST hTAG",
        "display": "NIST RM 8231 High-TAG",
        "title": "Ceramide Ring Trial — NIST RM 8231 High-TAG",
        "nativeId": "ILS_Ceramide_RingTrial_hTAG",
        "filename": "ils-ceramide-ring-trial-hTAG.json",
    },
    "DB": {
        "xlsx": "NIST DB",
        "display": "NIST RM 8231 Diabetic",
        "title": "Ceramide Ring Trial — NIST RM 8231 Diabetic",
        "nativeId": "ILS_Ceramide_RingTrial_DB",
        "filename": "ils-ceramide-ring-trial-DB.json",
    },
    "YAA": {
        "xlsx": "NIST YAA",
        "display": "NIST RM 8231 Young African-American",
        "title": "Ceramide Ring Trial — NIST RM 8231 Young African-American",
        "nativeId": "ILS_Ceramide_RingTrial_YAA",
        "filename": "ils-ceramide-ring-trial-YAA.json",
    },
}

# (Vendor, Instrument) exactly as written in data/definitions/instruments.csv,
# after .strip() on both fields.
INSTRUMENT_MODELS = {
    ("Agilent", "6460"): ("MS:1002445", "6460 Triple Quadrupole LC/MS"),
    ("Agilent", "6470"): ("MS:1003529", "6470 Triple Quadrupole LC/MS"),
    ("Agilent", "6490"): ("MS:1002446", "6490 Triple Quadrupole LC/MS"),
    ("Agilent", "6495"): ("MS:1003531", "6495 Triple Quadrupole LC/MS"),
    ("Bruker", "Bruker EVOQ Elite"): ("MS:1002297", "EVOQ Elite"),
    ("Sciex", "QTRAP 5500"): ("MS:1000931", "QTRAP 5500"),
    ("Sciex", "QTRAP 6500"): ("MS:1002581", "QTRAP 6500"),
    ("Sciex", "QTRAP 6500+"): ("MS:1002582", "QTRAP 6500+"),
    ("Sciex", "TripleTOF 6600"): ("MS:1002533", "TripleTOF 6600"),
    ("Thermo", "QE"): ("MS:1001911", "Q Exactive"),
    ("Thermo", "QE Plus"): ("MS:1002634", "Q Exactive Plus"),
    ("Thermo", "TSQ Quantiva"): ("MS:1002418", "TSQ Quantiva"),
    ("Waters", "Micromass Quattro Ultima"): ("MS:1000192", "Quattro Ultima"),
    ("Waters", "Synapt G2-Si"): ("MS:1002726", "SYNAPT G2-Si"),
    # Disambiguated from the original lab reports in data/reports/ (keyed by LabId):
    # report #32 names "Waters Xevo G2S QTOF", report #25 names "Xevo TQ-MS".
    ("Waters", "Xevo G2S"): ("MS:1002276", "Xevo G2-S QTof"),
    ("Waters", "Xevo TQ"): ("MS:1001790", "Xevo TQ MS"),
    ("Waters", "Xevo TQ-S"): ("MS:1001792", "Xevo TQ-S"),
}

VENDOR_MODELS = {
    "Sciex": ("MS:1000121", "SCIEX instrument model"),
    "Thermo": ("MS:1000483", "Thermo Fisher Scientific instrument model"),
    "Agilent": ("MS:1000490", "Agilent instrument model"),
    "Waters": ("MS:1000126", "Waters instrument model"),
    "Bruker": ("MS:1000122", "Bruker Daltonics instrument model"),
    "Shimadzu": ("MS:1000124", "Shimadzu instrument model"),
}

# LabIds that must NOT receive a model-level accession, mapped to the literal
# model string to place in `value`. Two lack a PSI-MS term; three are disputed
# between instruments.csv and the original lab reports. See design section 11.
INSTRUMENT_OVERRIDES = {
    "30": "Shimadzu LCMS-8060",     # PSI-MS has only LCMS-8060NX / -8060RX
    "35": "Sciex 4500 MD",          # PSI-MS has QTRAP 4500 / Triple Quad 4500, not the MD variant
    "38": "Waters Xevo TQ-S",       # report #38 says Xevo TQ-XS
    "29a": "Waters Synapt G2-Si",   # reports #29a/#29b appear swapped
    "29b": "Sciex QTRAP 6500",      # reports #29a/#29b appear swapped
}

ANALYZER_TYPES = {
    "QQQ": ("MS:1000081", "quadrupole"),
    "TOF": ("MS:1000084", "time-of-flight"),
    "Orbitrap": ("MS:1000484", "orbitrap"),
}

CHROMATOGRAPHY = {
    "RP": ("CHMO:0002302", "reversed-phase chromatography"),
    "FIA": ("MS:1000058", "flow injection analysis"),
    "SFC": ("MS:1003790", "supercritical fluid chromatography"),
}

PROTOCOLS = {
    "SOP": ("NCIT:C48443", "Standard Operating Procedure"),
    "OTHER": ("NCIT:C25294", "Laboratory Procedure"),
}

PUBLICATION_METADATA = {
    "title": ("Concordant inter-laboratory derived concentrations of ceramides in "
              "human plasma reference materials via authentic standards"),
    "authors": [
        "Federico Torta", "Nils Hoffmann", "Bo Burla", "Irina Alecu", "Makoto Arita",
        "Takeshi Bamba", "Steffany A. L. Bennett", "Justine Bertrand-Michel", "Britta Brügger",
        "Mónica P. Cala", "Dolores Camacho-Muñoz", "Antonio Checa", "Michael Chen",
        "Michaela Chocholoušková", "Michelle Cinel", "Emeline Chu-Van", "Benoit Colsch",
        "Cristina Coman", "Lisa Connell", "Bebiana C. Sousa", "Alex M. Dickens",
        "Maria Fedorova", "Finnur Freyr Eiríksson", "Hector Gallart-Ayala", "Mohan Ghorasaini",
        "Martin Giera", "Xue Li Guan", "Mark Haid", "Thomas Hankemeier", "Amy Harms",
        "Marcus Höring", "Michal Holčapek", "Thorsten Hornemann", "Chunxiu Hu",
        "Andreas J. Hülsmeier", "Kevin Huynh", "Christina M. Jones", "Julijana Ivanisevic",
        "Yoshihiro Izumi", "Harald C. Köfeler", "Sin Man Lam", "Mike Lange",
        "Jong Cheol Lee", "Gerhard Liebisch", "Katrice Lippa", "Andrea F. Lopez-Clavijo",
        "Malena Manzi", "Manuela R. Martinefski", "Raviswamy G. H. Math", "Satyajit Mayor",
        "Peter J. Meikle", "María Eugenia Monge", "Myeong Hee Moon", "Sneha Muralidharan",
        "Anna Nicolaou", "Thao Nguyen-Tran", "Valerie B. O'Donnell", "Matej Orešič",
        "Arvind Ramanathan", "Fabien Riols", "Daisuke Saigusa", "Tracey B. Schock",
        "Heidi Schwartz-Zimmermann", "Guanghou Shui", "Madhulika Singh",
        "Masatomo Takahashi", "Margrét Thorsteinsdóttir", "Noriyuki Tomiyasu",
        "Anthony Tournadre", "Hiroshi Tsugawa", "Victoria J. Tyrrell",
        "Grace van der Gugten", "Michael O. Wakelam", "Craig E. Wheelock",
        "Denise Wolrab", "Guowang Xu", "Tianrun Xu", "John A. Bowden",
        "Kim Ekroos", "Robert Ahrends", "Markus R. Wenk",
    ],
    "year": 2024,
    "doi": "10.1038/s41467-024-52087-x",
    "description": ("Nature Communications publication of the interlaboratory "
                    "ceramide ring trial study"),
}
