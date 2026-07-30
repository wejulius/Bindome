from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Set, Tuple, Optional
import re
import requests

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/{}.json"

CLASSES_14 = [
    "GPCRs",
    "Ion channels",
    "Transporters",
    "Intracellular signaling enzymes & adaptors",
    "Transcription & chromatin regulation",
    "DNA replication & chromosome biology",
    "RNA biology & translation",
    "Protein homeostasis & PTM control",
    "Metabolism & bioenergetics",
    "Cytoskeleton & motor proteins",
    "Membrane trafficking & organelles",
    "Cell junctions & adhesion",
    "Secreted & extracellular proteins",
    "Unassigned",
]   


@dataclass(frozen=True)
class ClassificationResult:
    uniprot_id: str
    class_name: str
    rule_id: str
    evidence: Dict[str, Any]


# -------------------------
# Fetch UniProt JSON
# -------------------------
def fetch_uniprot_entry(
    uniprot_id: str,
    timeout: int = 20,
    session: requests.Session | None = None,
) -> Dict[str, Any]:
    uid = uniprot_id.strip()
    client = session or requests
    r = client.get(
        UNIPROT_API.format(uid),
        headers={"Accept": "application/json"},
        timeout=timeout,
    )
    if r.status_code == 404:
        raise ValueError(f"UniProt ID not found: {uid}")
    r.raise_for_status()
    return r.json()


# -------------------------
# Structured extractors
# -------------------------
def _features(entry: Dict[str, Any]) -> List[Dict[str, Any]]:
    return entry.get("features", []) or []


def count_transmembranes(entry: Dict[str, Any]) -> int:
    return sum(1 for f in _features(entry) if (f.get("type") or "").lower() == "transmembrane")


def has_signal_peptide(entry: Dict[str, Any]) -> bool:
    return any((f.get("type") or "").lower() == "signal peptide" for f in _features(entry))


def has_gpi_anchor(entry: Dict[str, Any]) -> bool:
    return any((f.get("type") or "").lower() in {"gpi-anchor", "glycosylphosphatidylinositol-anchor"} for f in _features(entry))


def get_keywords(entry: Dict[str, Any]) -> Set[str]:
    out = set()
    for kw in entry.get("keywords", []) or []:
        name = kw.get("name")
        if isinstance(name, str) and name.strip():
            out.add(name.strip().lower())
    return out


def get_xrefs(entry: Dict[str, Any], db_name: str) -> List[Dict[str, Any]]:
    db_name_upper = db_name.upper()
    return [
        x for x in (entry.get("uniProtKBCrossReferences", []) or [])
        if (x.get("database") or "").upper() == db_name_upper
    ]


def get_go_terms(entry: Dict[str, Any]) -> Dict[str, Set[str]]:
    """
    UniProt GO xref properties often include:
      "F:protein kinase activity"
      "P:signal transduction"
      "C:nucleus"
    Returns dict: {"F": set(...), "P": set(...), "C": set(...)}
    """
    grouped = {"F": set(), "P": set(), "C": set()}
    for x in get_xrefs(entry, "GO"):
        for prop in x.get("properties", []) or []:
            v = prop.get("value")
            if not isinstance(v, str):
                continue
            v = v.strip()
            if len(v) >= 2 and v[1] == ":" and v[0] in grouped:
                grouped[v[0]].add(v[2:].strip().lower())
    return grouped


def get_interpro_ids(entry: Dict[str, Any]) -> Set[str]:
    ids = set()
    for x in get_xrefs(entry, "InterPro"):
        i = x.get("id")
        if isinstance(i, str) and i.strip():
            ids.add(i.strip().upper())
    return ids


def get_pfam_ids(entry: Dict[str, Any]) -> Set[str]:
    ids = set()
    for x in get_xrefs(entry, "Pfam"):
        i = x.get("id")
        if isinstance(i, str) and i.strip():
            ids.add(i.strip().upper())
    return ids


# -------------------------
# Protein name (high trust)
# -------------------------
def get_protein_name_text(entry: Dict[str, Any]) -> str:
    """
    Returns recommended + alternative names (NOT comments).
    This is high-trust per your preference.
    """
    texts: List[str] = []
    pd = entry.get("proteinDescription") or {}

    rec = pd.get("recommendedName") or {}
    if rec:
        fn = rec.get("fullName") or {}
        if isinstance(fn, dict) and isinstance(fn.get("value"), str):
            texts.append(fn["value"])

    for alt in (pd.get("alternativeNames") or []):
        fn = alt.get("fullName") or {}
        if isinstance(fn, dict) and isinstance(fn.get("value"), str):
            texts.append(fn["value"])

    # Some entries also have "submissionNames"
    for sub in (pd.get("submissionNames") or []):
        fn = sub.get("fullName") or {}
        if isinstance(fn, dict) and isinstance(fn.get("value"), str):
            texts.append(fn["value"])

    return " ".join(texts).lower()


# -------------------------
# Comments (lowest trust)
# -------------------------
def get_comments_text(entry: Dict[str, Any]) -> str:
    """
    Lowest trust. Only used as final tie-breaker before Unassigned.
    """
    parts: List[str] = []
    for c in entry.get("comments", []) or []:
        for t in c.get("texts", []) or []:
            v = t.get("value")
            if isinstance(v, str) and v.strip():
                parts.append(v)
    return " ".join(parts).lower()


def get_subcellular_localizations(entry: Dict[str, Any]) -> List[str]:
    locations: Set[str] = set()
    for c in entry.get("comments", []) or []:
        if (c.get("commentType") or "").upper() != "SUBCELLULAR LOCATION":
            continue
        for scl in c.get("subcellularLocations", []) or []:
            location = (scl.get("location") or {}).get("value")
            if isinstance(location, str) and location.strip():
                locations.add(location.strip())
    return sorted(locations)


# -------------------------
# Minimal, conservative regex patterns
# (applied to protein NAME first; comments last)
# -------------------------

RX_NAME = {
    # Structural / unambiguous
    "gpcr": re.compile(
        r"\b("
        r"g[- ]?protein[- ]?coupled receptor|gpcr|"
        r"adrenergic receptor|dopamine receptor|serotonin receptor|"
        r"histamine receptor|muscarinic receptor|opioid receptor|"
        r"chemokine receptor|melatonin receptor|prostaglandin receptor|"
        r"cannabinoid receptor|trace amine[- ]associated receptor|"
        r"frizzled receptor"
        r")\b",
        re.I,
    ),

    "ion_channel": re.compile(
        r"\b("
        r"ion channel|channel subunit|voltage[- ]gated|ligand[- ]gated|"
        r"potassium channel|sodium channel|calcium channel|chloride channel|"
        r"\bkcn[a-z0-9]+\b|\bscn\d+\b|\bcacn[a-z0-9]+\b|\bclcn\d+\b|"
        r"\btrp[cvma]\d*\b|"
        r"gap junction protein|connexin|pannexin"
        r")\b",
        re.I,
    ),

    "transporter": re.compile(
        r"\b("
        r"transporter|solute carrier|"
        r"\bslc\d+\b|"
        r"abc transporter|\babc[a-g]\d+\b|"
        r"pump|exchanger|symporter|antiporter|permease"
        r")\b",
        re.I,
    ),

    # Signaling enzymes & adaptors
    "kinase": re.compile(
        r"\b("
        r"kinase|mapk|map kinase|tyrosine kinase|serine/threonine kinase|"
        r"\bcdk\d+\b|\bjak\d+\b|\bpak\d+\b"
        r")\b",
        re.I,
    ),

    "phosphatase": re.compile(
        r"\b("
        r"phosphatase|"
        r"dual specificity phosphatase|"
        r"\bptp[a-z0-9]+\b"
        r")\b",
        re.I,
    ),

    "gtpase": re.compile(
        r"\b("
        r"gtpase|"
        r"\bras\b|\bran\b|"
        r"\brab\d+\b|"
        r"\barf\d+\b|"
        r"\brho[a-z0-9]*\b|\brac\d+\b|\bcdc42\b"
        r")\b",
        re.I,
    ),

    "cyclase": re.compile(r"\b(cyclase|adenylate cyclase|guanylate cyclase)\b", re.I),

    "phospholipase": re.compile(
        r"\b(phospholipase|phosphatidylinositol phospholipase|phospholipase a2|phospholipase c|phospholipase d)\b",
        re.I,
    ),

    # Transcription/chromatin
    "nuclear_receptor": re.compile(r"\b(nuclear receptor|steroid hormone receptor)\b", re.I),

    "transcription_factor": re.compile(
        r"\b("
        r"transcription factor|"
        r"homeobox|forkhead|"
        r"helix[- ]loop[- ]helix|"
        r"\bbzip\b|"
        r"ets domain|t[- ]box|"
        r"gata\b|sox\b|pax\b|"
        r"runx\b|tbx\b|"
        r"nf[- ]?kb"
        r")\b",
        re.I,
    ),

    "chromatin": re.compile(
        r"\b("
        r"chromatin|histone|"
        r"methyltransferase|demethylase|"
        r"acetyltransferase|deacetylase|"
        r"bromodomain|chromodomain|"
        r"polycomb|trithorax|"
        r"zinc finger protein|"
        r"\bzinc finger\b|"
        r"\bkrab\b|\bbtb\b"
        r")\b",
        re.I,
    ),

    # DNA/RNA
    "dna": re.compile(
        r"\b("
        r"dna polymerase|dna ligase|dna helicase|dna repair|"
        r"dna topoisomerase|"
        r"replication factor|"
        r"chromosome segregation"
        r")\b",
        re.I,
    ),

    "rna": re.compile(
        r"\b("
        r"ribosomal|translation|"
        r"spliceosome|splicing factor|"
        r"rna[- ]binding|rna binding|"
        r"ribonucleoprotein|"
        r"heterogeneous nuclear ribonucleoprotein|hnrnp"
        r")\b",
        re.I,
    ),

    # Proteostasis/PTM
    "ubiquitin": re.compile(
        r"\b("
        r"ubiquitin|e3 ligase|deubiquitin|"
        r"sumo|nedd8|isgyl|"
        r"ring finger protein|"
        r"tripartite motif|"
        r"\btrim\d+\b|"
        r"f[- ]box protein|"
        r"\bfbx[a-z0-9]+\b|"
        r"\bbtb domain\b"
        r")\b",
        re.I,
    ),

    "protease": re.compile(
        r"\b("
        r"protease|peptidase|"
        r"caspase|cathepsin|"
        r"metalloprotease|"
        r"serine protease|"
        r"ubiquitin[- ]specific protease|\busp\d+\b"
        r")\b",
        re.I,
    ),

    "chaperone": re.compile(
        r"\b("
        r"chaperone|heat shock|"
        r"\bhsp\d+\b|"
        r"chaperonin"
        r")\b",
        re.I,
    ),

    # Metabolism & bioenergetics
    "metabolic_enzyme": re.compile(
        r"\b("
        r"dehydrogenase|synthetase|synthase|isomerase|transferase|"
        r"oxidase|reductase|hydroxylase|mutase|lyase|hydrolase|"
        r"nadh dehydrogenase|"
        r"cytochrome c oxidase|cytochrome oxidase|"
        r"ubiquinol[- ]cytochrome c reductase|"
        r"succinate dehydrogenase|"
        r"atp synthase|"
        r"respiratory chain|ubiquinone"
        r")\b",
        re.I,
    ),

    # Structure/traffic/adhesion
    "cytoskeleton": re.compile(
        r"\b("
        r"actin|tubulin|microtubule|intermediate filament|"
        r"\bkeratin\b|"
        r"vimentin|desmin|nestin|peripherin|"
        r"neurofilament|lamin\b|"
        r"spectrin|filamin|"
        r"tropomyosin"
        r")\b",
        re.I,
    ),

    "motor": re.compile(r"\b(myosin|kinesin|dynein)\b", re.I),

    "trafficking": re.compile(
        r"\b("
        r"snare|clathrin|coatomer|copi|copii|"
        r"\brab\d+\b|"
        r"vesicle|"
        r"endosome|lysosome|golgi|"
        r"endoplasmic reticulum|"
        r"sorting nexin|exocyst"
        r")\b",
        re.I,
    ),

    "adhesion": re.compile(
        r"\b("
        r"integrin|cadherin|selectin|"
        r"tight junction|cell junction|"
        r"focal adhesion|"
        r"cell adhesion molecule|"
        r"\bicam\b|\bvcam\b|"
        r"immunoglobulin|ig[- ]like|"
        r"\bcd\d+\b|"
        r"t[- ]cell receptor|b[- ]cell receptor|"
        r"\bmhc\b|\bhla\b"
        r")\b",
        re.I,
    ),

    # Keep as a conservative NAME-level hint; in hierarchy gate with sp/gpi/extracellular GO/SCL.
    "secreted": re.compile(
        r"\b("
        r"cytokine|chemokine|interleukin|interferon|"
        r"growth factor|hormone|"
        r"defensin|cathelicidin|"
        r"complement|"
        r"collagen|laminin|fibronectin|proteoglycan|"
        r"secreted|extracellular|plasma protein"
        r")\b",
        re.I,
    ),
}

# Optional domain hints (InterPro/Pfam) — conservative examples.
# Add/extend as you wish.
INTERPRO_HINTS = {
    "GPCRs": {"IPR000276", "IPR017452"},  # 7TM GPCR, rhodopsin-like (examples)
    "Ion channels": {"IPR005821"},         # ion channel domain (example)
    "Transporters": {"IPR005828"},         # ABC transporter-like (example)
    "Transcription & chromatin regulation": {"IPR001356"},  # nuclear receptor (example)
    "Protein homeostasis & PTM control": {"IPR000608"},      # RING finger (example)
}

PFAM_HINTS = {
    "GPCRs": {"7TM_1", "7TM_2", "7TM_3"},
    "Ion channels": {"PF00520", "Kchannel", "Ion_trans"},
    "Transporters": {"ABC_tran", "PF00005"},
    "Transcription & chromatin regulation": {"Hormone_recep"},
    "Protein homeostasis & PTM control": {"RING"},
}


def classify_uniprot_id(
    uniprot_id: str,
    timeout: int = 20,
    session: requests.Session | None = None,
    entry: Dict[str, Any] | None = None,
) -> ClassificationResult:
    entry_data = entry if entry is not None else fetch_uniprot_entry(uniprot_id, timeout=timeout, session=session)

    kw = get_keywords(entry_data)                 # controlled vocabulary
    go = get_go_terms(entry_data)                 # structured GO text, grouped by F/P/C
    tm = count_transmembranes(entry_data)         # structured feature
    sp = has_signal_peptide(entry_data)           # structured feature
    gpi = has_gpi_anchor(entry_data)              # structured feature
    ipr = get_interpro_ids(entry_data)            # structured xref IDs
    pfam = get_pfam_ids(entry_data)               # structured xref IDs
    name = get_protein_name_text(entry_data)      # HIGH trust (per your request)
    comments = get_comments_text(entry_data)      # LOWEST trust
    scl = [s.lower() for s in get_subcellular_localizations(entry_data)]

    def kw_has(*needles: str) -> bool:
        return any(n.lower() in kw for n in needles)

    def go_has(aspect: str, *substrs: str) -> bool:
        hay = go.get(aspect, set())
        subs = [s.lower() for s in substrs]
        return any(any(ss in t for ss in subs) for t in hay)

    def ipr_hits(class_name: str) -> Set[str]:
        return set(INTERPRO_HINTS.get(class_name, set())).intersection(ipr)

    def pfam_hits(class_name: str) -> Set[str]:
        return set(PFAM_HINTS.get(class_name, set())).intersection(pfam)

    # Helper: match NAME first; only later consider COMMENTS
    def name_rx(key: str) -> bool:
        return bool(RX_NAME[key].search(name))

    def comments_rx(key: str) -> bool:
        return bool(RX_NAME[key].search(comments))

    def loc_has(*substrs: str) -> bool:
        needles = [s.lower() for s in substrs]
        return any(any(n in loc for n in needles) for loc in scl)

    # -------------------------
    # HIERARCHY: first match wins
    # -------------------------

    # 1) GPCRs (prefer topology + keyword/domain; then name)
    if (tm >= 7 and kw_has("g-protein coupled receptor")) or ipr_hits("GPCRs") or pfam_hits("GPCRs"):
        return ClassificationResult(uniprot_id, "GPCRs", "GPCR_STRUCT", {
            "tm_count": tm, "keyword_hit": "g-protein coupled receptor" if kw_has("g-protein coupled receptor") else None,
            "interpro_hits": sorted(ipr_hits("GPCRs")), "pfam_hits": sorted(pfam_hits("GPCRs"))
        })
    if tm >= 7 and name_rx("gpcr"):
        return ClassificationResult(uniprot_id, "GPCRs", "GPCR_NAME", {"tm_count": tm, "name_match": "gpcr"})

    # 2) Ion channels
    if kw_has("ion channel") or go_has("F", "ion channel activity", "channel activity") or ipr_hits("Ion channels") or pfam_hits("Ion channels"):
        return ClassificationResult(uniprot_id, "Ion channels", "ION_STRUCT", {
            "keyword_hit": "ion channel" if kw_has("ion channel") else None,
            "goF_channel_terms": sorted([t for t in go.get("F", set()) if "channel" in t]),
            "interpro_hits": sorted(ipr_hits("Ion channels")), "pfam_hits": sorted(pfam_hits("Ion channels"))
        })
    if name_rx("ion_channel"):
        return ClassificationResult(uniprot_id, "Ion channels", "ION_NAME", {"name_match": "ion_channel"})

    # 3) Transporters
    if kw_has("transporter", "solute carrier", "atp-binding cassette transporter") or go_has("F", "transporter activity") or ipr_hits("Transporters") or pfam_hits("Transporters"):
        return ClassificationResult(uniprot_id, "Transporters", "TRANSPORT_STRUCT", {
            "keyword_hits": sorted(list(kw.intersection({"transporter", "solute carrier", "atp-binding cassette transporter"}))),
            "goF_transporter_terms": sorted([t for t in go.get("F", set()) if "transporter" in t]),
            "interpro_hits": sorted(ipr_hits("Transporters")), "pfam_hits": sorted(pfam_hits("Transporters"))
        })
    if name_rx("transporter"):
        return ClassificationResult(uniprot_id, "Transporters", "TRANSPORT_NAME", {"name_match": "transporter"})

    # 4) Intracellular signaling enzymes & adaptors (GO MF / keywords first; then name)
    if (
        go_has("F", "kinase activity", "protein kinase activity", "phosphatase activity", "gtpase activity", "guanyl-nucleotide exchange factor activity")
        or kw_has("kinase", "phosphatase", "small gtpase")
    ):
        return ClassificationResult(uniprot_id, "Intracellular signaling enzymes & adaptors", "SIGNAL_STRUCT", {
            "goF_hits": sorted([t for t in go.get("F", set()) if any(k in t for k in ["kinase", "phosphatase", "gtpase", "exchange factor"])]),
            "keyword_hits": sorted(list(kw.intersection({"kinase", "phosphatase", "small gtpase"}))),
        })
    if any(name_rx(k) for k in ["kinase", "phosphatase", "gtpase", "cyclase", "phospholipase"]):
        return ClassificationResult(uniprot_id, "Intracellular signaling enzymes & adaptors", "SIGNAL_NAME", {
            "name_matches": [k for k in ["kinase", "phosphatase", "gtpase", "cyclase", "phospholipase"] if name_rx(k)]
        })

    # 5) Transcription & chromatin regulation (GO MF/keywords/domains first; then name)
    if (
        go_has("F", "dna-binding transcription factor activity", "transcription factor activity", "chromatin binding")
        or kw_has("transcription", "dna-binding", "chromatin regulator", "histone modification")
        or ipr_hits("Transcription & chromatin regulation")
        or pfam_hits("Transcription & chromatin regulation")
    ):
        return ClassificationResult(uniprot_id, "Transcription & chromatin regulation", "TX_STRUCT", {
            "goF_hits": sorted([t for t in go.get("F", set()) if any(k in t for k in ["transcription", "chromatin", "dna-binding"])]),
            "keyword_hits": sorted(list(kw.intersection({"transcription", "dna-binding", "chromatin regulator", "histone modification"}))),
            "interpro_hits": sorted(ipr_hits("Transcription & chromatin regulation")),
            "pfam_hits": sorted(pfam_hits("Transcription & chromatin regulation")),
        })
    if any(name_rx(k) for k in ["nuclear_receptor", "transcription_factor", "chromatin"]):
        return ClassificationResult(uniprot_id, "Transcription & chromatin regulation", "TX_NAME", {
            "name_matches": [k for k in ["nuclear_receptor", "transcription_factor", "chromatin"] if name_rx(k)]
        })

    # 6) DNA replication & chromosome biology
    if (
        go_has("P", "dna replication", "dna repair", "chromosome segregation")
        or go_has("F", "dna polymerase activity", "dna helicase activity", "dna ligase activity")
    ):
        return ClassificationResult(uniprot_id, "DNA replication & chromosome biology", "DNA_STRUCT", {
            "goP_hits": sorted([t for t in go.get("P", set()) if "dna" in t or "chromosome" in t]),
            "goF_hits": sorted([t for t in go.get("F", set()) if t.startswith("dna") or "dna " in t]),
        })
    if name_rx("dna"):
        return ClassificationResult(uniprot_id, "DNA replication & chromosome biology", "DNA_NAME", {"name_match": "dna"})

    # 7) RNA biology & translation
    if (
        go_has("P", "mRNA splicing", "translation", "ribosome biogenesis")
        or go_has("F", "structural constituent of ribosome", "rna binding")
        or kw_has("ribosomal protein", "rna-binding", "translation")
    ):
        return ClassificationResult(uniprot_id, "RNA biology & translation", "RNA_STRUCT", {
            "goP_hits": sorted([t for t in go.get("P", set()) if any(k in t for k in ["splicing", "translation", "ribosome"])]),
            "goF_hits": sorted([t for t in go.get("F", set()) if any(k in t for k in ["ribosome", "rna binding"])]),
            "keyword_hits": sorted(list(kw.intersection({"ribosomal protein", "rna-binding", "translation"}))),
        })
    if name_rx("rna"):
        return ClassificationResult(uniprot_id, "RNA biology & translation", "RNA_NAME", {"name_match": "rna"})

    # 8) Protein homeostasis & PTM control
    if (
        kw_has("ubiquitin", "ubiquitin-like conjugation", "protease", "proteasome", "chaperone", "autophagy")
        or go_has("P", "protein ubiquitination", "proteolysis", "protein folding")
        or go_has("F", "ubiquitin-protein transferase activity", "peptidase activity")
        or ipr_hits("Protein homeostasis & PTM control")
        or pfam_hits("Protein homeostasis & PTM control")
    ):
        return ClassificationResult(uniprot_id, "Protein homeostasis & PTM control", "HOMEOSTASIS_STRUCT", {
            "keyword_hits": sorted(list(kw.intersection({
                "ubiquitin", "ubiquitin-like conjugation", "protease", "proteasome", "chaperone", "autophagy"
            }))),
            "goP_hits": sorted([t for t in go.get("P", set()) if any(k in t for k in ["ubiquitin", "proteolysis", "folding"])]),
            "goF_hits": sorted([t for t in go.get("F", set()) if any(k in t for k in ["ubiquitin", "peptidase"])]),
            "interpro_hits": sorted(ipr_hits("Protein homeostasis & PTM control")),
            "pfam_hits": sorted(pfam_hits("Protein homeostasis & PTM control")),
        })
    if any(name_rx(k) for k in ["ubiquitin", "protease", "chaperone"]):
        return ClassificationResult(uniprot_id, "Protein homeostasis & PTM control", "HOMEOSTASIS_NAME", {
            "name_matches": [k for k in ["ubiquitin", "protease", "chaperone"] if name_rx(k)]
        })

    # 9) Metabolism & bioenergetics
    if go_has("P", "metabolic process", "biosynthetic process", "catabolic process") and go_has(
        "F",
        "oxidoreductase activity",
        "transferase activity",
        "ligase activity",
        "isomerase activity",
        "lyase activity",
        "hydrolase activity",
    ):
        return ClassificationResult(uniprot_id, "Metabolism & bioenergetics", "METABOLISM_GO", {
            "goP_hits": sorted([t for t in go.get("P", set()) if "metabolic" in t or "biosynthetic" in t or "catabolic" in t]),
            "goF_enzyme_hits": sorted([t for t in go.get("F", set()) if any(k in t for k in ["oxidoreductase", "transferase", "ligase", "isomerase", "lyase", "hydrolase"])]),
        })
    if name_rx("metabolic_enzyme") and go_has("P", "metabolic process"):
        # require metabolic GO P to avoid random "transferase" names from non-metabolic contexts
        return ClassificationResult(uniprot_id, "Metabolism & bioenergetics", "METABOLISM_NAME+GO", {
            "name_match": "metabolic_enzyme",
            "goP_hits": sorted([t for t in go.get("P", set()) if "metabolic" in t]),
        })

    # 10) Cytoskeleton & motor proteins
    if go_has("C", "cytoskeleton") or go_has("P", "cytoskeleton organization"):
        return ClassificationResult(uniprot_id, "Cytoskeleton & motor proteins", "CYTO_GO", {
            "goC_hits": sorted([t for t in go.get("C", set()) if "cytoskeleton" in t]),
            "goP_hits": sorted([t for t in go.get("P", set()) if "cytoskeleton" in t]),
        })
    if any(name_rx(k) for k in ["cytoskeleton", "motor"]):
        return ClassificationResult(uniprot_id, "Cytoskeleton & motor proteins", "CYTO_NAME", {
            "name_matches": [k for k in ["cytoskeleton", "motor"] if name_rx(k)]
        })

    # 11) Membrane trafficking & organelles
    if go_has("P", "vesicle-mediated transport", "endocytosis") or go_has("C", "golgi apparatus", "endosome", "lysosome", "endoplasmic reticulum"):
        return ClassificationResult(uniprot_id, "Membrane trafficking & organelles", "TRAFFICKING_GO", {
            "goP_hits": sorted([t for t in go.get("P", set()) if any(k in t for k in ["vesicle", "endocytosis"])]),
            "goC_hits": sorted([t for t in go.get("C", set()) if any(k in t for k in ["golgi", "endosome", "lysosome", "endoplasmic reticulum"])]),
        })
    if name_rx("trafficking"):
        return ClassificationResult(uniprot_id, "Membrane trafficking & organelles", "TRAFFICKING_NAME", {"name_match": "trafficking"})

    # 12) Cell junctions & adhesion
    if go_has("P", "cell adhesion") or go_has("C", "cell junction"):
        return ClassificationResult(uniprot_id, "Cell junctions & adhesion", "ADHESION_GO", {
            "goP_hits": sorted([t for t in go.get("P", set()) if "adhesion" in t]),
            "goC_hits": sorted([t for t in go.get("C", set()) if "junction" in t]),
        })
    if name_rx("adhesion"):
        return ClassificationResult(uniprot_id, "Cell junctions & adhesion", "ADHESION_NAME", {"name_match": "adhesion"})

    # 13) Secreted & extracellular proteins
    # Prefer structured secretion evidence: signal peptide/GPI and extracellular GO C/keywords.
    if (sp or gpi) and (go_has("C", "extracellular region", "extracellular space") or kw_has("secreted")):
        return ClassificationResult(uniprot_id, "Secreted & extracellular proteins", "SECRETED_STRUCT", {
            "signal_peptide": sp, "gpi_anchor": gpi,
            "goC_hits": sorted([t for t in go.get("C", set()) if "extracellular" in t]),
            "keyword_hit": "secreted" if kw_has("secreted") else None,
            "tm_count": tm,
        })
    # Protein name as fallback (still above comments)
    if (sp or gpi or tm == 0) and name_rx("secreted"):
        return ClassificationResult(uniprot_id, "Secreted & extracellular proteins", "SECRETED_NAME", {
            "signal_peptide": sp, "gpi_anchor": gpi, "tm_count": tm, "name_match": "secreted"
        })

    # 14) Subcellular localization (weak fallback)
    # Applied only after stronger structured/name rules above.
    if loc_has("extracellular", "cell surface", "secreted", "extracellular matrix"):
        return ClassificationResult(uniprot_id, "Secreted & extracellular proteins", "SECRETED_LOC_WEAK", {
            "subcellular_localizations": scl,
        })
    if loc_has("cell junction", "desmosome", "tight junction", "adherens junction", "focal adhesion"):
        return ClassificationResult(uniprot_id, "Cell junctions & adhesion", "ADHESION_LOC_WEAK", {
            "subcellular_localizations": scl,
        })
    if loc_has("golgi", "endosome", "lysosome", "endoplasmic reticulum", "vesicle", "secretory granule"):
        return ClassificationResult(uniprot_id, "Membrane trafficking & organelles", "TRAFFICKING_LOC_WEAK", {
            "subcellular_localizations": scl,
        })
    if loc_has("cytoskeleton", "microtubule", "actin", "myofibril", "sarcomere", "centrosome", "spindle"):
        return ClassificationResult(uniprot_id, "Cytoskeleton & motor proteins", "CYTO_LOC_WEAK", {
            "subcellular_localizations": scl,
        })
    if loc_has("nucleolus", "spliceosome", "ribosome", "ribonucleoprotein"):
        return ClassificationResult(uniprot_id, "RNA biology & translation", "RNA_LOC_WEAK", {
            "subcellular_localizations": scl,
        })
    if loc_has("chromosome", "replication fork"):
        return ClassificationResult(uniprot_id, "DNA replication & chromosome biology", "DNA_LOC_WEAK", {
            "subcellular_localizations": scl,
        })
    if loc_has("nucleus", "chromatin"):
        return ClassificationResult(uniprot_id, "Transcription & chromatin regulation", "TX_LOC_WEAK", {
            "subcellular_localizations": scl,
        })
    if loc_has("mitochond", "peroxisome"):
        return ClassificationResult(uniprot_id, "Metabolism & bioenergetics", "METABOLISM_LOC_WEAK", {
            "subcellular_localizations": scl,
        })

    # 15) COMMENTS (lowest priority) — only if you explicitly want a last-ditch rescue
    # Keep this extremely conservative: only accept when comments contain a direct controlled-like phrase AND nothing else matched.
    if (sp or gpi) and re.search(r"\bsecreted\b|\bextracellular\b", comments, flags=re.I):
        return ClassificationResult(uniprot_id, "Secreted & extracellular proteins", "SECRETED_COMMENTS_LAST", {
            "signal_peptide": sp, "gpi_anchor": gpi, "comment_hit": True
        })

    return ClassificationResult(uniprot_id, "Unassigned", "UNASSIGNED", {
        "tm_count": tm, "signal_peptide": sp, "gpi_anchor": gpi,
        "keywords": sorted(list(kw))[:50],  # truncate for logging
        "goF_count": len(go.get("F", set())), "goP_count": len(go.get("P", set())), "goC_count": len(go.get("C", set())),
        "interpro_count": len(ipr), "pfam_count": len(pfam),
    })