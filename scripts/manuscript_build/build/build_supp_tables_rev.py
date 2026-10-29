"""Build Additional file 2 (Tables S1-S18) for the revised Human Genomics submission from results/ and results/revision/.

Each sheet: table title in row 1, then one or more blocks (block subtitle, bold header, data). Header row of the
first block is frozen. P-like columns are stored as numbers rounded to 3 significant digits and shown as 0.00E+00.
Run on the Windows host:  python build_supp_tables_rev.py
"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import re, math, json
from pathlib import Path
import numpy as np, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

R = Path(f"{PROJ_ROOT}/results")
OUT = Path(f"{PROJ_ROOT}/submission/revision/additional_files/Additional_file_2_Tables_S1-S18.xlsx")
LOG = []

P_COL = re.compile(r"(^p$|^p[_.]|_p$|_p_|^P$|^P_|_P$|pval|p_value|fdr|bonf|^pC|_pC|binom_p|^Q_p$|steiger_p|presso_global_p|p_HEIDI|p_SMR|p_eQTL|p_GWAS|min_p|p_dec|p_smr|p_amd|p_csvd|^p_trend$|p_interaction|PP\.H|PP4|joint_PP|lead_pip|coloc_snp_PP|SNP\.PP)",
                   re.I)
SCI_COL = re.compile(r"(^p$|^p[_.]|_p$|_p_|^P$|^P_|_P$|pval|p_value|fdr|bonf|^pC|_pC|binom_p|^Q_p$|steiger_p|presso_global_p|p_HEIDI|p_SMR|p_eQTL|p_GWAS|min_p|p_dec|p_smr|p_amd|p_csvd|^p_trend$|p_interaction)",
                     re.I)


def rd(f, **kw):
    f = R / f
    d = pd.read_csv(f, sep="\t", **kw)
    d.columns = [str(c).strip() for c in d.columns]
    if d.columns[0] in ("", "Unnamed: 0"):
        d = d.rename(columns={d.columns[0]: "gene" if "gene" not in d.columns else "item"})
    return d


def sig(x, n):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))) or x == 0:
        return x
    return float(f"{x:.{n}g}")


def tidy(d):
    d = d.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            n = 3 if SCI_COL.search(c) else 4
            d[c] = d[c].map(lambda v: sig(v, n) if pd.notna(v) else v)
    return d


class Book:
    def __init__(self):
        self.wb = Workbook(); self.wb.remove(self.wb.active); self.contents = []

    def sheet(self, num, title, desc, blocks, coldefs="", source=""):
        ws = self.wb.create_sheet(f"Table S{num}")
        ws.cell(1, 1, f"Table S{num}. {title}").font = Font(bold=True, size=12)
        r = 3; first_header = None; nrows = 0
        for sub, df in blocks:
            df = tidy(df)
            ws.cell(r, 1, sub).font = Font(bold=True, italic=True); r += 1
            for j, c in enumerate(df.columns, 1):
                cell = ws.cell(r, j, c); cell.font = Font(bold=True)
                cell.fill = PatternFill("solid", fgColor="E7E6E6"); cell.alignment = Alignment(wrap_text=True, vertical="top")
            if first_header is None: first_header = r
            sci = [j for j, c in enumerate(df.columns, 1) if SCI_COL.search(c)]
            for row in df.itertuples(index=False):
                r += 1
                for j, v in enumerate(row, 1):
                    try:
                        if v is pd.NA or v is pd.NaT or (isinstance(v, float) and math.isnan(v)): v = None
                    except TypeError:
                        pass
                    if isinstance(v, (np.integer,)): v = int(v)
                    if isinstance(v, (np.floating,)): v = float(v)
                    if isinstance(v, (np.bool_,)): v = bool(v)
                    if isinstance(v, str) and len(v) > 32000: v = v[:32000] + "...(truncated)"
                    cell = ws.cell(r, j, v)
                    if j in sci and isinstance(v, float) and v != 0: cell.number_format = "0.00E+00"
            nrows += len(df); r += 3
            LOG.append(f"S{num} | {sub} | {len(df)} rows x {df.shape[1]} cols")
        for j in range(1, ws.max_column + 1):
            ws.column_dimensions[get_column_letter(j)].width = 16
        ws.column_dimensions["A"].width = 24
        if first_header: ws.freeze_panes = ws.cell(first_header + 1, 1)
        self.contents.append(dict(table=f"Table S{num}", title=title, description=desc, n_rows=nrows,
                                  blocks="; ".join(b[0] for b in blocks), column_definitions=coldefs, source_files=source))

    def save(self, abbrev):
        ws = self.wb.create_sheet("Contents", 0)
        ws.cell(1, 1, "Additional file 2. Supplementary Tables S1-S18").font = Font(bold=True, size=13)
        ws.cell(2, 1, "Genetic overlap between age-related macular degeneration and stroke / cerebral small vessel disease. "
                      "All effect estimates are on the log-odds (binary traits) or standardised beta (quantitative traits) scale unless stated; "
                      "P values are two-sided and shown to three significant digits.")
        hdr = ["Table", "Title", "Description", "Rows", "Blocks (sub-tables)", "Column definitions", "Source (results/)"]
        for j, h in enumerate(hdr, 1):
            c = ws.cell(4, j, h); c.font = Font(bold=True); c.fill = PatternFill("solid", fgColor="E7E6E6")
        for i, x in enumerate(self.contents, 5):
            for j, k in enumerate(["table", "title", "description", "n_rows", "blocks", "column_definitions", "source_files"], 1):
                c = ws.cell(i, j, x[k]); c.alignment = Alignment(wrap_text=True, vertical="top")
        r = 5 + len(self.contents) + 2
        ws.cell(r, 1, "Abbreviations").font = Font(bold=True)
        for k, v in abbrev:
            r += 1; ws.cell(r, 1, k); ws.cell(r, 2, v)
        for col, w in zip("ABCDEFG", [11, 42, 60, 8, 40, 70, 40]):
            ws.column_dimensions[col].width = w
        ws.freeze_panes = "A5"
        OUT.parent.mkdir(parents=True, exist_ok=True)
        self.wb.save(OUT)


RV = R / "revision"


def rr(f, **kw):
    return rd("revision/" + f, **kw)


def r2(f, **kw):
    return rd("revision2/" + f, **kw)


def tag(d, **cols):
    d = d.copy()
    for i, (k, v) in enumerate(cols.items()):
        d.insert(i, k, v)
    return d


B = Book()

# ======================= S1 data sources
src_cols = ["role", "dataset", "phenotype_definition", "ancestry", "n_cases", "n_controls", "n_total", "genome_build", "accession_or_url", "PMID"]
S1 = pd.DataFrame([
    ["AMD (primary exposure)", "FinnGen R12 H7_AMD", "Any AMD: ICD-10 H35.30 or H35.31 (ICD-9 3625A/3625B), hospital discharge or cause-of-death register, onset age >= 50 (see endpoint block)", "Finnish", 12495, 461686, 474181, "GRCh38", "FinnGen R12 public results (finngen_R12_H7_AMD)", "36653562"],
    ["AMD (primary exposure)", "FinnGen R12 DRY_AMD", "Dry AMD: ICD-10 H35.30 / ICD-9 3625A (non-exudative AMD of any registered stage)", "Finnish", 8570, 329258, 337828, "GRCh38", "finngen_R12_DRY_AMD", "36653562"],
    ["AMD (primary exposure)", "FinnGen R12 WET_AMD", "Wet AMD: ICD-10 H35.31 / ICD-9 3625B", "Finnish", 6699, 331070, 337769, "GRCh38", "finngen_R12_WET_AMD", "36653562"],
    ["AMD (replication)", "IAMDGC (Fritsche 2016)", "Advanced AMD (geographic atrophy and/or choroidal neovascularisation); effect sizes approximated from Z; P < 1e-300 recovered from mantissa/exponent", "Mainly European", 16144, 17832, 33976, "GRCh38 (harmonised)", "GCST003219", "26691988"],
    ["AMD (instruments, sensitivity)", "IAMDGC 52 conditionally independent variants", "Conditionally independent AMD variants (scheme B_52)", "Mainly European", None, None, None, "GRCh38", "GWAS Catalog associations for GCST003219", "26691988"],
    ["AMD (exploratory)", "Early AMD (Winkler 2020)", "Early / intermediate AMD (includes IAMDGC and UK Biobank participants)", "European", 14034, 91214, 105248, "GRCh38 (harmonised)", "GCST010723", "32843070"],
    ["Stroke (outcome/exposure)", "GIGASTROKE EUR AS", "Any stroke", "European", 73652, 1234808, 1308460, "GRCh38 (harmonised)", "GCST90104539", "36180795"],
    ["Stroke (outcome/exposure)", "GIGASTROKE EUR IS", "Ischaemic stroke", "European", 62100, 1234808, 1296908, "GRCh38 (harmonised)", "GCST90104540", "36180795"],
    ["Stroke (outcome/exposure)", "GIGASTROKE EUR CES", "Cardioembolic stroke (TOAST)", "European", 10804, 1234808, 1245612, "GRCh38 (harmonised)", "GCST90104541", "36180795"],
    ["Stroke (outcome/exposure)", "GIGASTROKE EUR LAS", "Large artery stroke (TOAST)", "European", 6399, 1234808, 1241207, "GRCh38 (harmonised)", "GCST90104542", "36180795"],
    ["Stroke (outcome/exposure)", "GIGASTROKE EUR SVS", "Small vessel stroke (TOAST)", "European", 6811, 1234808, 1241619, "GRCh38 (harmonised)", "GCST90104543", "36180795"],
    ["Stroke (FinnGen-free replication)", "MEGASTROKE EUR AS", "Any stroke (no FinnGen/FINRISK; other Finnish cohorts contributed)", "European", 40585, 406111, 446696, "GRCh38 (harmonised)", "GCST006906", "29531354"],
    ["Stroke (FinnGen-free replication)", "MEGASTROKE EUR LAS", "Large artery stroke", "European", 4373, 297290, 301663, "GRCh38 (harmonised)", "GCST006907", "29531354"],
    ["Stroke (FinnGen-free replication)", "MEGASTROKE EUR IS", "Ischaemic stroke", "European", 34217, 406111, 440328, "GRCh38 (harmonised)", "GCST006908", "29531354"],
    ["Stroke (FinnGen-free replication)", "MEGASTROKE EUR SVS", "Small vessel stroke", "European", 5386, 343560, 348946, "GRCh38 (harmonised)", "GCST006909", "29531354"],
    ["Stroke (FinnGen-free replication)", "MEGASTROKE EUR CES", "Cardioembolic stroke", "European", 7193, 355468, 362661, "GRCh38 (harmonised)", "GCST006910", "29531354"],
    ["Stroke (cross-ancestry)", "GIGASTROKE EAS AS", "Any stroke", "East Asian", 27413, 237242, 264655, "GRCh37 (rsID matched)", "GCST90104544", "36180795"],
    ["Stroke (cross-ancestry)", "GIGASTROKE EAS IS", "Ischaemic stroke", "East Asian", 19032, 237242, 256274, "GRCh37 (rsID matched)", "GCST90104545", "36180795"],
    ["Stroke (cross-ancestry)", "GIGASTROKE EAS LAS", "Large artery stroke", "East Asian", 1735, 237242, 238977, "GRCh37 (rsID matched)", "GCST90104547", "36180795"],
    ["Stroke (cross-ancestry)", "GIGASTROKE EAS SVS", "Small vessel stroke", "East Asian", 5532, 237242, 242774, "GRCh37 (rsID matched)", "GCST90104548", "36180795"],
    ["CSVD", "UK Biobank BIG40 IDP 1437", "White matter hyperintensity volume (MRI)", "European (UK Biobank)", None, None, 33224, "GRCh37", "https://open.win.ox.ac.uk/ukbiobank/big40/ (phenotype 1437)", "33875891"],
    ["CSVD", "Lacunar stroke (Traylor 2021)", "MRI-confirmed or clinically defined lacunar stroke", "European", 6030, 219389, 225419, "GRCh38 (harmonised)", "GCST90014122", "33773637"],
    ["CSVD", "EPVS white matter (Duperron 2023)", "Extensive perivascular spaces, white matter", "European", 9317, 29281, 38598, "GRCh38 (harmonised)", "GCST90244151", "37069360"],
    ["CSVD", "EPVS basal ganglia (Duperron 2023)", "Extensive perivascular spaces, basal ganglia", "European", 8950, 29953, 38903, "GRCh38 (harmonised)", "GCST90244153", "37069360"],
    ["CSVD", "EPVS hippocampus (Duperron 2023)", "Extensive perivascular spaces, hippocampus", "European", 9163, 29708, 38871, "GRCh38 (harmonised)", "GCST90244155", "37069360"],
    ["CSVD", "WM-MD (Persyn 2020)", "White matter mean diffusivity, first principal component", "European (UK Biobank)", None, None, 17467, "GRCh37", "GCST010103", "32358547"],
    ["CSVD", "WM-FA (Persyn 2020)", "White matter fractional anisotropy, first principal component", "European (UK Biobank)", None, None, 17663, "GRCh37", "GCST010102", "32358547"],
    ["Haemorrhagic stroke", "ICH, UK Biobank WGS", "Intracerebral haemorrhage (ICD-10 I61), non-Finnish European; no overlap with the AMD exposures", "European", 2407, 456033, 458440, "GRCh38 (harmonised)", "GCST90473584", "40770095"],
    ["Locus-level (3p14)", "cIMT (Franceschini 2018)", "Common carotid intima-media thickness", "European", None, None, 71128, "GRCh37/38", "GCST007436", "30510157"],
    ["Locus-level (3p14)", "Carotid plaque (Franceschini 2018)", "Presence of carotid plaque", "European", 21540, 26894, 48434, "GRCh37/38", "GCST007435", "30510157"],
    ["Locus-level (3p14)", "CAD (Aragam 2022)", "Coronary artery disease", "Mainly European", 181522, 984168, 1165690, "GRCh38 (harmonised)", "GCST90132314", "36474045"],
    ["Locus-level / MVMR", "Type 2 diabetes (Xue 2018)", "Type 2 diabetes (UK Biobank + DIAGRAM)", "European", 62892, 596424, 659316, "GRCh38 (harmonised)", "GCST006867", "30054458"],
    ["Locus-level / MVMR", "WHRadjBMI (Pulit 2019)", "Waist-to-hip ratio adjusted for BMI", "European", None, None, 694649, "GRCh38 (harmonised)", "GCST008999", "30239722"],
    ["Locus-level / MVMR", "SBP (Evangelou 2018)", "Systolic blood pressure", "European", None, None, 757601, "GRCh37 (position matched)", "GCST006624", "30224653"],
    ["Locus-level / MVMR", "DBP (Evangelou 2018)", "Diastolic blood pressure", "European", None, None, 757601, "GRCh37 (position matched)", "GCST006630", "30224653"],
    ["Locus-level (3p14)", "LDL cholesterol (GLGC 2021)", "LDL cholesterol, European", "European", None, None, 1320016, "GRCh38 (harmonised)", "GCST90239658", "34887591"],
    ["Locus-level (3p14)", "HDL cholesterol (GLGC 2021)", "HDL cholesterol, European", "European", None, None, 1320016, "GRCh38 (harmonised)", "GCST90239652", "34887591"],
    ["Locus-level (FinnGen endpoints)", "FinnGen R12 T2D", "Type 2 diabetes, definitions combined", "Finnish", 82878, 403489, 486367, "GRCh38", "FinnGen R12 (remote tabix)", "36653562"],
    ["Locus-level (FinnGen endpoints)", "FinnGen R12 KELA_DIAB_INSUL_EXMORE", "Diabetes, insulin treatment (Kela reimbursement)", "Finnish", 77557, 403479, 481036, "GRCh38", "FinnGen R12", "36653562"],
    ["Locus-level (FinnGen endpoints)", "FinnGen R12 I9_STR", "Stroke, excluding subarachnoid haemorrhage", "Finnish", 34110, 450023, 484133, "GRCh38", "FinnGen R12", "36653562"],
    ["Locus-level (FinnGen endpoints)", "FinnGen R12 C_STROKE", "Any stroke", "Finnish", 53492, 360342, 413834, "GRCh38", "FinnGen R12", "36653562"],
    ["Locus-level (FinnGen endpoints)", "FinnGen R12 I9_ICH", "Intracerebral haemorrhage", "Finnish", 5112, 450016, 455128, "GRCh38", "FinnGen R12", "36653562"],
    ["Locus-level (FinnGen endpoints)", "FinnGen R12 I9_CHD", "Major coronary heart disease event", "Finnish", 56650, 443698, 500348, "GRCh38", "FinnGen R12", "36653562"],
    ["Locus-level (FinnGen endpoints)", "FinnGen R12 I9_CORATHER", "Coronary atherosclerosis", "Finnish", 63307, 416171, 479478, "GRCh38", "FinnGen R12", "36653562"],
    ["Locus-level (FinnGen endpoints)", "FinnGen R12 I9_HYPTENS", "Hypertension", "Finnish", 154630, 345634, 500264, "GRCh38", "FinnGen R12", "36653562"],
    ["Fine-mapping", "FinnGen R12 SuSiE fine-mapping", "In-sample SuSiE credible sets for H7_AMD, DRY_AMD, WET_AMD, T2D and insulin-treated diabetes", "Finnish", None, None, None, "GRCh38", "FinnGen R12 public fine-mapping release", "36653562"],
    ["LD reference", "FinnGen SISu v4.2 LD", "In-sample LD (r2, D') for key 10q26 variant pairs", "Finnish", None, None, None, "GRCh38", "FinnGen LD server (SISu v4.2)", "36653562"],
    ["LD reference", "1000 Genomes Phase 3 EUR / EAS", "LD reference panels (EUR n = 503; EAS for East Asian conditional models)", "European; East Asian", None, None, 503, "GRCh37/38", "https://www.internationalgenome.org/", "26432245"],
    ["cis-eQTL", "eQTLGen phase I", "Whole-blood cis-eQTL (BESD for SMR)", "Mainly European", None, None, 31684, "GRCh37", "https://www.eqtlgen.org/", "34475573"],
    ["cis-eQTL", "GTEx v8 (via eQTL Catalogue)", "49 tissues incl. aorta, coronary and tibial artery, brain, whole blood, fibroblasts; aortic SuSiE credible sets (QTD000131, n = 387)", "Mainly European", None, None, None, "GRCh38", "eQTL Catalogue (remote tabix)", "32913098; 34493866"],
    ["cis-eQTL", "ROSMAP DLPFC (via eQTL Catalogue)", "Dorsolateral prefrontal cortex", "European", None, None, 560, "GRCh38", "eQTL Catalogue", "28869584; 34493866"],
    ["cis-eQTL", "BrainSeq DLPFC (via eQTL Catalogue)", "Dorsolateral prefrontal cortex", "European/African American", None, None, 479, "GRCh38", "eQTL Catalogue", "31174959; 34493866"],
    ["cis-eQTL", "EyeGEx", "Neural retina", "European", None, None, 406, "GRCh37", "Zenodo record 18331341", "30742112"],
    ["cis-pQTL", "UKB-PPP (Sun 2023), discovery", "Olink Explore; one primary cis-pQTL per protein (1,792 proteins)", "European (UK Biobank)", None, None, None, "GRCh38", "Nature 2023 Supplementary Table 9", "37794186"],
    ["cis-pQTL", "deCODE (Ferkingstad 2021), replication", "SomaScan; cis-pQTL rank 1 (1,449 proteins); HTRA1 pQTL", "Icelandic", None, None, None, "GRCh38", "Nat Genet 2021 Supplementary Table 2", "34857953"],
    ["Single-cell", "RPE-choroid atlas (CELLxGENE)", "scRNA-seq of human RPE and choroid (182,151 cells; normal donors)", "Mixed", None, None, 182151, "-", "CELLxGENE dataset 5ed90bc0-12f3-4354-b314-0ce0dbb96458; doi:10.64898/2026.07.23.740352", "-"],
    ["Single-cell", "Human retina, periphery (Cowan 2020)", "scRNA-seq of peripheral retina (34,723 cells)", "Mixed", None, None, 34723, "-", "CELLxGENE dataset ccfdcad0-7104-46b9-addf-fd66a2a15907", "32946783"],
    ["Single-cell", "Human Brain Cell Atlas v1.0, vascular supercluster (Siletti 2023)", "snRNA-seq of brain vascular cells (9,932 nuclei)", "Mixed", None, None, 9932, "-", "CELLxGENE dataset 0d7f4c06-a6bd-47d2-a42b-3a7196704f77", "37824663"],
    ["Single-cell", "Brain vascular multi-omics (Reid 2025)", "snRNA-seq of brain vascular cells (65,479 nuclei)", "Mixed", None, None, 65479, "-", "CELLxGENE dataset 203025fe-fa99-4d57-81da-458ed8f0c334", "40730185"],
    ["Observational", "NHANES 2005-2008", "Retinal photograph-graded AMD (age >= 40) and self-reported stroke", "US population", 249, 4787, 5036, "-", "https://wwwn.cdc.gov/nchs/nhanes/", "21220632"],
], columns=src_cols)
for c in ("n_cases", "n_controls", "n_total"):
    S1[c] = S1[c].astype("Int64")
FG_END = pd.DataFrame([
    ["H7_AMD", "Age-related macular degeneration (dry or wet)", "H35.30 or H35.31", "3625A or 3625B", "Hospital discharge; cause of death", "Event age >= 50 (age < 50 '_EARLY' endpoints excluded)", "Diabetic retinopathy, retinal detachment, other macular disease (MACULA_NAS), all glaucoma endpoints", 12495, 461686, "Union of DRY_AMD and WET_AMD"],
    ["DRY_AMD", "Dry (non-exudative) AMD", "H35.30", "3625A", "Hospital discharge; cause of death", "Event age >= 50", "As H7_AMD; controls also require death age >= 50 if deceased", 8570, 329258, "Covers non-exudative AMD of any registered stage, not only geographic atrophy"],
    ["WET_AMD", "Wet (exudative) AMD", "H35.31", "3625B", "Hospital discharge; cause of death", "Event age >= 50", "As DRY_AMD", 6699, 331070, "About 2,774 cases carry both dry and wet codes (8,570 + 6,699 - 12,495)"],
], columns=["endpoint", "label", "ICD-10", "ICD-9", "registers", "age_condition", "control_exclusions", "n_cases", "n_controls", "note"])
B.sheet(1, "Data sources and FinnGen AMD endpoint definitions",
        "GWAS, fine-mapping, QTL, LD-reference, single-cell and survey data sources with phenotype definitions, sample sizes, accession numbers and PMIDs; FinnGen R12 AMD endpoint definitions (Risteys R12).",
        [("Data sources", S1), ("FinnGen R12 AMD endpoint definitions", FG_END)],
        "role = use in the study; n_cases/n_controls for binary traits, n_total for quantitative traits or QTL datasets; PMID = PubMed ID of the source publication (verified). FinnGen stroke endpoint counts (I9_STR, C_STROKE, I9_ICH) derived from N x case fraction. Sample overlap: FinnGen contributed <= 205,672 participants (<= 15.7%) to GIGASTROKE EUR (Table S5).",
        "scripts/30_prepare_gwas.R; 45_csvd_gwas_prep.R; rev_mr_00_prep.R; 76_rev_extra.R; revision/adamts9/traits_used.tsv; revision/adamts9/finngen_phewas_rs4611812.tsv; revision/extra/SUMMARY.md")

# ======================= S2 NHANES (unchanged)
B.sheet(2, "NHANES 2005-2008: sample selection, baseline characteristics and survey-weighted associations between AMD and stroke",
        "Participant flow; characteristics by stroke status; weighted stroke prevalence by AMD category; models M1-M4; sensitivity and subgroup analyses.",
        [("Participant flow", rd("nhanes/flowchart.tsv")),
         ("Baseline characteristics by stroke status (survey-weighted)", rd("nhanes/table1_by_stroke.tsv")),
         ("Weighted stroke prevalence by AMD category", rd("nhanes/stroke_prevalence_by_amd.tsv")),
         ("Main models (M1 unadjusted; M2 + demographics; M3 + socioeconomic and lifestyle; M4 + vascular risk factors)", rd("nhanes/main_models.tsv")),
         ("Sensitivity analyses", rd("nhanes/sensitivity.tsv")),
         ("Subgroup analyses with multiplicative interaction", rd("nhanes/subgroups.tsv"))],
        "OR/lo/hi = odds ratio and 95% CI from survey-weighted logistic regression; p_trend = trend across no/early/late AMD; p_interaction = Wald test for the product term. Model 4 n = 4,843 (complete full-model covariates).",
        "nhanes/")

# ======================= S3 instruments
rev_inst = []
for f in sorted((R / "mr").glob("instruments_*.tsv")):
    if any(k in f.stem for k in ("AMD",)):
        continue
    x = pd.read_csv(f, sep="\t"); x.insert(0, "exposure", f.stem.replace("instruments_", "")); rev_inst.append(x)
for f in sorted((R / "csvd").glob("instruments_*.tsv")):
    if any(k in f.stem for k in ("AMD", "SVS", "WMH")):
        continue
    x = pd.read_csv(f, sep="\t"); x.insert(0, "exposure", f.stem.replace("instruments_", "")); rev_inst.append(x)
rev_inst = pd.concat(rev_inst, ignore_index=True)
B.sheet(3, "Genetic instruments for Mendelian randomization",
        "AMD instruments under three schemes (A_dist: one lead per +/- 1 Mb, primary; B_52: 52 conditionally independent IAMDGC variants; C_orig: LD clumping r2 < 0.001 / 10 Mb), instruments per locus, outcome attrition and proxies, early-AMD instruments and instruments for reverse-direction MR.",
        [("AMD instruments, all schemes", rr("mr/instruments_all_schemes.tsv")),
         ("Instruments per locus (10q26, ADAMTS9, APOE, C3, CFH, lipid loci, other)", rr("mr/instruments_per_locus.tsv")),
         ("Instrument attrition per outcome (direct, proxied, lost)", rr("mr/instrument_attrition.tsv")),
         ("Outcome proxy summary (union of instruments)", rr("mr/outcome_proxy_summary.tsv")),
         ("Early AMD instruments (Winkler 2020)", rr("extra/instruments_EARLYAMD_W20.tsv")),
         ("Instruments for reverse MR (stroke, WMH and CSVD exposures; LD clumping)", rev_inst)],
        "scheme A_dist = lead variants (P < 5e-8, F > 10) of distance-defined loci, one per +/- 1 Mb window (primary); B_52 = correlated-instrument IVW with signed 1000G LD; C_orig = clumping r2 < 0.001 within 10 Mb. Proxies: 1000G r2 > 0.8 within +/- 500 kb with phased allele mapping. n_lost = instruments with neither the variant nor a proxy in the outcome GWAS.",
        "revision/mr/instruments_all_schemes.tsv; instruments_per_locus.tsv; instrument_attrition.tsv; outcome_proxy_summary.tsv; revision/extra/instruments_EARLYAMD_W20.tsv; mr/instruments_*.tsv; csvd/instruments_*.tsv")

# ======================= S4 MR results
old = rd("mr/mr_all_methods.tsv")
rev_dir = old[~old["exposure"].str.contains("AMD") & (old["analysis"] == "main")].drop(columns=["id.exposure", "id.outcome"], errors="ignore")
cs = rd("csvd/mr_all_methods.tsv").drop(columns=["id.exposure", "id.outcome"], errors="ignore")
cs = cs[~cs["outcome"].isin(["SVS", "WMH"]) | ~cs["exposure"].str.contains("AMD")]
B.sheet(4, "Mendelian randomization estimates (all directions, instrument schemes and methods)",
        "Forward MR of four AMD phenotypes on GIGASTROKE and MEGASTROKE stroke outcomes and WMH for three instrument schemes, ten locus-exclusion sets and all estimators (IVW, MR-Egger, weighted median, weighted mode, MR-RAPS, correlated IVW/Egger), with power and equivalence; reverse MR; CSVD extension traits; ICH and early AMD.",
        [("Forward MR, AMD -> stroke (GIGASTROKE, MEGASTROKE) and WMH", rr("mr/mr_revised_all.tsv")),
         ("Reverse MR, stroke subtypes and WMH -> AMD (LD-clumped instruments; main analysis)", rev_dir),
         ("MR with CSVD extension traits (lacunar stroke, EPVS, white-matter MD/FA), both directions (LD-clumped instruments)", cs),
         ("MR with ICH and early AMD (IVW with FDR)", rr("extra/mr_ich_earlyamd_ivw.tsv").drop(columns=["id.exposure", "id.outcome"], errors="ignore")),
         ("MR with ICH and early AMD (all methods)", rr("extra/mr_ich_earlyamd_all_methods.tsv").drop(columns=["id.exposure", "id.outcome"], errors="ignore"))],
        "b/se per unit log-odds of exposure liability (binary outcomes: OR = exp(b); WMH in SD); exclusion = locus set removed (+/- 2 Mb): no_10q26, no_ADAMTS9, no_CFH, no_C3, no_APOE, no_major3 (CFH, 10q26, C3), no_major3_APOE, no_major3_lipid, no_all; mde_logOR / mde_OR_up / mde_OR_down = minimum detectable effect at 80% power; tost_p = two one-sided tests of equivalence (bounds OR 0.95-1.05; +/- 0.05 SD for WMH); equivalent = tost_p < 0.05. MEGA_* = MEGASTROKE EUR outcomes.",
        "revision/mr/mr_revised_all.tsv; mr/mr_all_methods.tsv; csvd/mr_all_methods.tsv; revision/extra/mr_ich_earlyamd_*.tsv")

# ======================= S5 MR sensitivity
mra = rr("mr/mr_revised_all.tsv")
pw = mra[(mra["method"] == "IVW") & (mra["scheme"] == "A_dist") & (mra["exclusion"].isin(["main", "no_10q26", "no_major3", "no_all"]))]
pw = pw[["exposure", "scheme", "outcome", "exclusion", "OR", "lo", "hi", "p", "nsnp", "mde_logOR", "mde_OR_down", "mde_OR_up", "tost_p", "equivalent"]]
d1 = rd("mr/mr_diagnostics.tsv"); d1.insert(0, "analysis_set", "AMD and stroke/WMH (LD-clumped instruments)")
d2 = rd("csvd/mr_diagnostics.tsv"); d2.insert(0, "analysis_set", "CSVD extension (LD-clumped instruments)")
d3 = rr("extra/mr_ich_earlyamd_diagnostics.tsv"); d3.insert(0, "analysis_set", "ICH and early AMD")
B.sheet(5, "Mendelian randomization sensitivity analyses: locus attribution, clustering, pleiotropy, sample overlap, power and equivalence",
        "Leave-one-locus-out estimates (+/- 2 Mb), locus contributions to the IVW estimate, MR-Clust, CAUSE, MRlap correction for sample overlap, quantified sample overlap, power and equivalence tests (primary on estimates without 10q26), translation of the equivalence bounds to the scale of manifest AMD, and conventional diagnostics (Cochran Q, MR-Egger intercept, MR-PRESSO, Steiger).",
        [("Leave-one-locus-out IVW estimates (wide)", rr("mr/mr_leave_one_locus_out_wide.tsv")),
         ("Locus contributions to IVW estimates", rr("mr/mr_locus_contributions.tsv")),
         ("MR-Clust cluster summary", rr("mr/mrclust_summary.tsv")), ("MR-Clust variant assignments", rr("mr/mrclust_assignments.tsv")),
         ("CAUSE (FinnGen AMD -> IS, SVS)", rr("mr/cause_results.tsv")),
         ("MRlap: sample-overlap-corrected estimates", rr("mr/mrlap_overlap_corrected.tsv")),
         ("Sample overlap quantification", rr("mr/sample_overlap_quantification.tsv")),
         ("Power (minimum detectable OR at 80%) and equivalence (TOST), primary IVW", pw),
         ("Equivalence by exclusion set (no_10q26 primary; no_major3 / no_all sensitivity)", r2("equivalence_tost_by_exclusion.tsv")),
         ("Translation of the equivalence bounds to the causal OR of manifest AMD (Burgess and Labrecque)", r2("equivalence_bounds.tsv")),
         ("Positive controls: AMD instruments in the IAMDGC; systolic blood pressure on stroke and CAD", rd("story/positive_controls.tsv")),
         ("Per-locus Wald ratios of AMD loci on stroke, WMH and CAD (pathway annotation)", rd("story/per_locus_effects.tsv")),
         ("MR partitioned by AMD pathway (complement, 10q26, lipid, other) and test against the causal expectation", rd("story/pathway_mr.tsv")),
         ("AMD liability and coronary artery disease (Aragam 2022)", rd("story/amd_cad_mr.tsv")),
         ("Scaling of exposure liability to early-AMD liability (kappa)", rd("revision5/kappa_earlyamd_scale.tsv")),
         ("Re-specified comparison with the effect expected under causality (kappa scaling, exact derivative, propagated uncertainty)", rd("revision5/causal_expectation_respecified.tsv")),
         ("Early-AMD liability as exposure", rd("revision5/earlyamd_exposure_mr.tsv")),
         ("MR by AMD definition (registry, advanced, early) and instrument scheme", rd("revision6/phenotype_matched_mr.tsv")),
         ("Difference between the expected (scaled observational) and observed MR effects", rd("revision6/benchmark_difference_ci.tsv")),
         ("Sign-concordance: binomial versus block-permutation inference", rd("revision6/sign_concordance_inference.tsv")),
         ("Positive control: atrial fibrillation liability on stroke", rd("revision5/positive_control_af.tsv")),
         ("Complement: per-variant Wald ratios (functional variants and locus leads)", rd("revision5/complement_per_variant.tsv")),
         ("Complement: heterogeneity within the group and gene-based membership", rd("revision5/complement_heterogeneity.tsv")),
         ("Complement: gene annotation of instruments", rd("revision5/complement_instrument_genes.tsv")),
         ("Complement: conditionally independent IAMDGC variants, correlated-instrument IVW", rd("revision5/complement_iamdgc_cor_ivw.tsv")),
         ("Earlier specification of the benchmark test (unscaled; superseded)", rd("revision4/causal_expectation_test.tsv")),
         ("Bound-free equivalence: smallest symmetric bound at which equivalence can be declared", rd("revision4/equivalence_boundfree.tsv")),
         ("TOST P value as a function of the equivalence bound", rd("revision4/tost_curve.tsv")),
         ("Conventional diagnostics", pd.concat([d1, d2, d3], ignore_index=True)),
         ("Steiger directionality on the liability scale (LD-clumped instruments)", rd("mr/mr_steiger_liability.tsv"))],
        "b_main / b_no_* = IVW estimate (log-odds) with all instruments or after removing the locus (+/- 2 Mb); weight_share = share of IVW weight; cluster_class = null / positive / negative MR-Clust cluster; p_causal_vs_sharing = CAUSE test of causal vs sharing-only model; corrected_b = MRlap estimate; max_overlap_fraction = upper bound of shared participants.",
        "revision/mr/mr_leave_one_locus_out_wide.tsv; mr_locus_contributions.tsv; mrclust_*.tsv; cause_results.tsv; mrlap_overlap_corrected.tsv; sample_overlap_quantification.tsv; mr_revised_all.tsv; mr/mr_diagnostics.tsv; csvd/mr_diagnostics.tsv; revision/extra/mr_ich_earlyamd_diagnostics.tsv; mr/mr_steiger_liability.tsv")

# ======================= S6 LDSC
fa = rr("mr/fdr_all_layers.tsv")
ld_fdr = fa[fa["layer"].str.startswith("LDSC")]
B.sheet(6, "Genome-wide genetic correlation (LD score regression)",
        "SNP heritability and cross-trait genetic correlations between AMD and stroke, CSVD, ICH and early AMD; FDR for the primary LDSC family (FinnGen AMD, 36 tests); cross-trait intercept matrix used for LAVA.",
        [("FDR, primary family (FinnGen AMD x 9 traits)", ld_fdr),
         ("SNP heritability, AMD and stroke/WMH", rd("ldsc/ldsc_h2.tsv")), ("Genetic correlations, AMD and stroke/WMH", rd("ldsc/ldsc_rg_all.tsv")),
         ("SNP heritability, CSVD extension", rd("csvd/ldsc_h2.tsv")), ("Genetic correlations, CSVD extension", rd("csvd/ldsc_rg.tsv")),
         ("SNP heritability, ICH and early AMD", rr("extra/ldsc_h2_ich_earlyamd.tsv")), ("Genetic correlations, ICH and early AMD", rr("extra/ldsc_rg_ich_earlyamd.tsv")),
         ("Cross-trait LDSC intercept matrix", pd.read_csv(R / "csvd/ldsc_intercept_matrix.tsv", sep="\t"))],
        "rg = genetic correlation; gcov_int = cross-trait intercept; fdr = BH within the LDSC layer. IAMDGC LDSC estimates are invalid (Z-derived effect sizes) and were not interpreted.",
        "ldsc/; csvd/ldsc_*.tsv; revision/extra/ldsc_*_ich_earlyamd.tsv; revision/mr/fdr_all_layers.tsv")

# ======================= S7 coloc
B.sheet(7, "Pairwise GWAS-GWAS colocalisation and prior sensitivity",
        "coloc.abf for AMD and stroke lead loci (+/- 500 kb), annotated shared loci, prior sensitivity over p12 = 1e-5 to 1e-6 at 10q26 for all AMD and stroke datasets, and prior dependence of the earlier results with the fraction of stroke evidence on AMD credible-set variants.",
        [("Shared loci (PP.H4 > 0.5), annotated", rd("coloc/coloc_shared_loci_annotated.tsv")),
         ("10q26: coloc.abf prior sensitivity (p12)", rr("10q26/coloc_abf_p12_sensitivity.tsv")),
         ("Prior dependence of earlier colocalisation results", rr("10q26/earlier_coloc_prior_dependence.tsv")),
         ("All locus x trait-pair tests (default priors)", rd("coloc/coloc_gwas_gwas.tsv"))],
        "PP.H3 = distinct causal variants; PP.H4 = shared causal variant (single-causal-variant assumption; default p1 = p2 = 1e-4, p12 = 1e-5); ratio = PP.H4/PP.H3; f_shared = fraction of the partner trait's posterior evidence on variants of the AMD 95% credible set.",
        "coloc/; revision/10q26/coloc_abf_p12_sensitivity.tsv; revision/10q26/earlier_coloc_prior_dependence.tsv")

# ======================= S8 LAVA
B.sheet(8, "Local genetic correlation (LAVA)",
        "All 448 bivariate local genetic correlations with within-trait and global FDR, summary by unique LD block, 10q26 sensitivity analyses (pairwise inputs, MEGASTROKE outcomes, sample-overlap scenarios, merged blocks), tests per trait pair, and follow-up of the concordant blocks (lead-variant look-ups and coloc.abf).",
        [("Bivariate local rg with within-trait and global FDR", rr("mr/lava_bivariate_with_global_fdr.tsv")),
         ("Summary by unique LD block", rr("mr/lava_unique_blocks.tsv")),
         ("10q26 sensitivity analyses", rr("mr/lava_10q26_sensitivity.tsv")),
         ("Tests per trait pair", rd("csvd/lava_counts.tsv")),
         ("Concordant blocks: lead-variant look-ups and colocalisation", r2("concordant_blocks_followup.tsv")),
         ("Concordant blocks: lead variants", r2("concordant_blocks_leads.tsv"))],
        "rho = local genetic correlation (95% CI rho.lower-rho.upper); p_fdr = BH within each stroke/CSVD trait; fdr_global = BH over all 448 tests; scenario = pairwise input / MEGASTROKE / overlap assumption (0, LDSC intercept, 0.1, 0.2) / merged blocks; for WM-FA higher values indicate better integrity. Univariate local heritability for all 2,495 blocks: results/csvd/lava_univariate.tsv.",
        "revision/mr/lava_bivariate_with_global_fdr.tsv; lava_unique_blocks.tsv; lava_10q26_sensitivity.tsv; csvd/lava_counts.tsv")

# ======================= S9 concordance
B.sheet(9, "Effect-direction concordance of AMD lead variants in stroke and CSVD traits",
        "One lead variant per distance-defined AMD locus aligned to its risk allele; concordant fraction, binomial test and block sign-flip permutation (20,000 permutations; leads within 5 Mb flipped together) with FDR; variant-level data.",
        [("Concordance tests", rr("mr/concordance_distance_pruned.tsv")), ("Variant-level data", rr("mr/concordance_variant_level.tsv"))],
        "n_loci = AMD loci available for the trait; frac_same = concordant fraction; perm_p_frac / perm_p_sumz = permutation P for the fraction / sum of signed z; fdr_perm_sumz = BH over traits; z_amd, z_o = signed z for AMD and the outcome; block = permutation block.",
        "revision/mr/concordance_*.tsv")

# ======================= S10 key variants
B.sheet(10, "Effects of key shared variants across AMD, stroke (European and East Asian), CSVD, ICH, early AMD and cardiometabolic traits",
        "Variant-level effects of 10q26 (rs3750848, rs10490924, rs36212733, rs11200638, rs2672592, rs72631113), 3p14 (rs4611812, rs6795735, rs17676309), APOE (rs429358, rs7412), C4 (rs1269852, rs1270942), ABO (rs947073006), SH2B3 (rs3184504), CFH and C3 variants.",
        [("10q26 variants in all datasets (AMD-risk allele; corrected IAMDGC)", rr("10q26/key_variants_all_datasets.tsv")),
         ("Key variants in ICH, early AMD and East Asian stroke", rr("extra/key_variants_ich_earlyamd_eas.tsv")),
         ("Key variants in CSVD traits (FinnGen AMD-risk allele)", rd("csvd/key_variants_AMD_FG.tsv")),
         ("Key variants in CSVD traits (IAMDGC AMD-risk allele)", rd("csvd/key_variants_AMD_IAMDGC.tsv")),
         ("Key variants in AMD and GIGASTROKE EUR", rd("tissue_eqtl/key_variant_gwas_lookup.tsv")),
         ("3p14 rs4611812-C effects across traits", rr("adamts9/rs4611812C_effects_all_traits.tsv")),
         ("10q26 variants and EPVS (primary sample-size meta-analysis Z, AMD-risk allele)", r2("epvs_corrected.tsv")),
         ("Summary of loci shared by AMD and stroke or CSVD outside 10q26 (moved from main Table 3)", rd("revision7/other_shared_loci_summary.tsv"))],
        "risk / risk_allele / amd_risk_allele / allele = effect allele (AMD-risk allele in FinnGen unless stated); beta_risk, b_risk, b_csvd, b_allele = log-odds or SD; OR_risk = exp(b_risk); Neff_rel = effective sample size relative to the maximum in the region; rs4611812-C = the AMD-protective, aortic ADAMTS9-raising allele. IAMDGC rows in the GIGASTROKE look-up predate the IAMDGC P-value correction; use the first block for IAMDGC at 10q26.",
        "revision/10q26/key_variants_all_datasets.tsv; revision/extra/key_variants_ich_earlyamd_eas.tsv; csvd/key_variants_*.tsv; tissue_eqtl/key_variant_gwas_lookup.tsv; revision/adamts9/rs4611812C_effects_all_traits.tsv")

# ======================= S11 conditional
tq = rd("tissue_eqtl/conditional_two_snp.tsv")
q12 = tq[tq["region"].astype(str).str.contains("12q24")]
cj = rd("tissue_eqtl/cojo_cond.tsv")
B.sheet(11, "Conditional analyses at 10q26, 12q24 and 6p21 (C4)",
        "10q26: the same conditional models for AMD and stroke (two-variant models in both directions, joint three-variant models, 2-df test of rs2672592/rs72631113 given rs3750848) with ancestry-matched LD (FinnGen SISu v4.2, 1000G EUR, 1000G EAS; r +/- 0.03 and 0.1) in GIGASTROKE, MEGASTROKE, the approximate GIGASTROKE remainder and FinnGen; projected proportionality test (European data only); GCTA-COJO; 12q24: rs3184504 vs SH2B3 eQTL lead rs587914; 6p21: C4 eQTL variants conditioned on CFB/C2 AMD leads.",
        [("10q26 two-variant conditional models, ancestry-matched LD", rr("10q26/conditional_two_snp_matched_LD.tsv")),
         ("10q26 joint three-variant models (rs3750848, rs2672592, rs72631113)", rr("10q26/joint_three_snp_model.tsv")),
         ("10q26 symmetric conditional models for AMD and stroke (2-df side | lead, 1-df lead | side, joint models)", r2("10q26_symmetric_conditional.tsv")),
         ("10q26 LD sensitivity of the conditional tests", r2("10q26_r_sensitivity.tsv")),
         ("10q26 conditional models in GIGASTROKE minus MEGASTROKE and FinnGen (approximate)", r2("10q26_giga_remainder_conditional.tsv")),
         ("10q26 Pan-UK Biobank: symmetric conditional models with UK Biobank in-sample LD", rd("revision3/panukb/panukb_symmetric_conditional.tsv")),
         ("10q26 Pan-UK Biobank: LD sensitivity", rd("revision3/panukb/panukb_r_sensitivity.tsv")),
         ("10q26 Pan-UK Biobank: proportionality tests", rd("revision3/panukb/panukb_proportionality.tsv")),
         ("10q26 Pan-UK Biobank: four-variant joint models (with coronary signal rs4237541)", rd("revision3/panukb/fourvar_joint_models.tsv")),
         ("10q26 proportionality test (AMD joint effects projected through stroke LD)", (lambda d: d[~d["stroke"].astype(str).str.contains("EAS")])(r2("10q26_proportionality_projected.tsv"))),
         ("GCTA-COJO stepwise selection (1000G EUR; not interpretable for AMD)", rd("tissue_eqtl/cojo_slct.tsv")),
         ("GCTA-COJO conditional analysis (10q26 and 12q24)", cj),
         ("12q24 two-variant conditional models (1000G EUR)", q12),
         ("6p21 C4: lead variants", rr("extra/c4_leads.tsv")),
         ("6p21 C4: conditional on CFB/C2 AMD leads", rr("extra/c4_conditional.tsv"))],
        "r = LD correlation from ld_source; r_shift = sensitivity shift (-0.1, 0, +0.1); pC_b_given_a = conditional P of SNP b given SNP a; z_joint / p_joint = joint-model estimates; p_side_given_lead_2df = joint 2-df test of rs2672592 and rs72631113 given rs3750848; p_lead_given_side_1df = rs3750848 given both markers; chi2_projected / p_projected = proportionality test with the FinnGen AMD joint effects projected through the stroke LD (ridge = regularisation of the LD matrix; East Asian stroke data not used); the earlier unprojected proportionality results are superseded; pC_c4_given_cond = conditional P of the C4 eQTL variant.",
        "revision/10q26/conditional_two_snp_matched_LD.tsv; joint_three_snp_model.tsv; proportionality_test.tsv; tissue_eqtl/cojo_*.tsv; tissue_eqtl/conditional_two_snp.tsv; revision/extra/c4_*.tsv")

# ======================= S12 fine-mapping
fgcs = pd.read_csv(RV / "10q26/finngen_susie_credsets_10q26.tsv", sep="\t", header=None,
                   names=["trait", "region", "cs", "cs_log10bf", "cs_avg_r2", "cs_min_r2", "low_purity", "cs_size"])
ldj = json.loads((RV / "10q26/ld_finngen_sisu42_pairs.json").read_text())
ldf = pd.DataFrame([{"snp1": k.split("|")[0], "snp2": k.split("|")[1], "r2": v.get("r2"), "dprime": v.get("dprime")} for k, v in ldj.items()])
B.sheet(12, "Fine-mapping at 10q26: FinnGen in-sample SuSiE, external-panel SuSiE-RSS, coloc.bf_bf, LD structure and effective sample size",
        "FinnGen R12 in-sample SuSiE credible sets and PIPs of key variants for any, dry and wet AMD; FinnGen SISu v4.2 in-sample LD; SuSiE-RSS with 1000G EUR and kriging_rss diagnostics; fine-mapping-based colocalisation (coloc.bf_bf) of FinnGen AMD credible sets with stroke; LD (r, D') in 1000G EUR and Ensembl EUR/FIN/EAS; regional effective N for GIGASTROKE; FinnGen stroke endpoints, heterogeneity between stroke datasets, FinnGen power and winner's-curse correction.",
        [("FinnGen in-sample SuSiE credible sets (10q26 region)", fgcs),
         ("FinnGen in-sample SuSiE: key variants", rr("10q26/finngen_susie_key_variants.tsv")),
         ("FinnGen SISu v4.2 in-sample LD for key variant pairs", ldf),
         ("SuSiE-RSS with 1000G EUR LD: credible sets", rd("tissue_eqtl/susie_10q26_credible_sets.tsv")),
         ("SuSiE-RSS: LD-mismatch diagnostics (kriging_rss)", rd("tissue_eqtl/susie_10q26_ld_mismatch_diagnostics.tsv")),
         ("coloc.bf_bf: FinnGen AMD SuSiE vs stroke", rr("10q26/coloc_bf_bf_AMDsusie_vs_stroke.tsv")),
         ("LD, 1000G EUR (r and D')", rr("10q26/ld_1000G_EUR_r_dprime.tsv")),
         ("LD, Ensembl EUR / FIN / EAS", rr("10q26/ld_ensembl_EUR_FIN_EAS.tsv")),
         ("Regional effective sample size (GIGASTROKE, FinnGen stroke)", rr("10q26/effective_N_region.tsv")),
         ("FinnGen stroke endpoints at 10q26 (I9_STR, C_STROKE, I9_TIA)", r2("finngen_extra_endpoints_10q26.tsv")),
         ("FinnGen vs other cohorts: heterogeneity, power and winner's curse", r2("finngen_heterogeneity_power.tsv")),
         ("Heterogeneity across MEGASTROKE, GIGASTROKE remainder and FinnGen", r2("10q26_stroke_heterogeneity_3parts.tsv")),
         ("Pan-UK Biobank phenotypes used at 10q26", rd("revision3/panukb/panukb_phenotypes.tsv")),
         ("Pan-UK Biobank: key variants (European ancestry)", rd("revision3/panukb/panukb_10q26_variants.tsv")),
         ("Pan-UK Biobank: UK Biobank in-sample LD for key variant pairs", rd("revision3/panukb/panukb_ld_keypairs.tsv")),
         ("Pan-UK Biobank: SuSiE with in-sample LD (summary)", rd("revision3/panukb/panukb_susie_summary.tsv")),
         ("Pan-UK Biobank: SuSiE credible-set members", rd("revision3/panukb/panukb_susie_cs_members.tsv")),
         ("Heterogeneity across MEGASTROKE, GIGASTROKE remainder, FinnGen and UK Biobank", rd("revision3/panukb/heterogeneity_with_other_datasets.tsv")),
         ("10q26 haplotype-by-trait matrix of joint effects (AMD, CAD, IS, SVS, WMH)", rd("story/haplotype_trait_matrix.tsv")),
         ("Sample-overlap audit: cohort overlap matrix", rd("revision6/overlap_matrix.tsv")),
         ("Sample-overlap audit: null-variant Z-statistic correlations between datasets", rd("revision6/null_z_correlations.tsv")),
         ("GIGASTROKE remainder: per-variant construction by inverse-variance subtraction (information shares, before and after)", rd("revision7/remainder_construction.tsv")),
         ("GIGASTROKE remainder: two-haplotype joint model before subtraction, for each subtracted dataset and after subtraction", rd("revision7/remainder_joint_before_after.tsv")),
         ("10q26 stroke: meta-analysis modeling sample overlap (GLS) and strictly non-overlapping subsets", rd("revision6/10q26_stroke_gls_overlap.tsv")),
         ("10q26 stroke: datasets with in-sample LD only (FinnGen, UK Biobank)", rd("revision6/insample_only_pool.tsv")),
         ("10q26 stroke: dataset-level heterogeneity table", rd("revision6/heterogeneity_main_table.tsv")),
         ("10q26 coronary disease: independent datasets (FinnGen, CARDIoGRAMplusC4D 2015)", rd("revision6/cad_independent_validation.tsv")),
         ("10q26 coronary disease: pooled independent datasets (with and without UK Biobank)", rd("revision6/cad_independent_pooled.tsv")),
         ("10q26: power of the conditional and difference tests", rd("revision6/contrast_power.tsv")),
         ("10q26: simulation of shared versus distinct causal haplotypes", rd("revision6/simulation_10q26.tsv")),
         ("10q26 stroke: haplotype contrast and Hartung-Knapp/REML random-effects inference", rd("revision5/10q26_stroke_contrast_hksj.tsv")),
         ("10q26 stroke: dataset-level joint effects and contrasts", rd("revision5/10q26_stroke_contrast_datasets.tsv")),
         ("Random-effects meta-analysis of joint 10q26 effects across four European stroke datasets", rd("revision4/10q26_stroke_random_effects.tsv")),
         ("Random-effects meta-analysis: leave-one-dataset-out", rd("revision4/10q26_stroke_loo.tsv")),
         ("Random-effects meta-analysis: dataset-level joint effects", rd("revision4/10q26_stroke_datasets.tsv"))],
        "cs = credible set index (-1 = not in a credible set); PIP = posterior inclusion probability; cs_log10bf, cs_avg_r2, cs_min_r2, low_purity, cs_size = FinnGen credible-set summary columns (header assigned from the FinnGen release format); hit1/hit2 = credible-set lead variants; Neff = 4/(1/ncase + 1/ncontrol) per variant. 1000G SuSiE-RSS fits for AMD are invalid (LD mismatch) and shown for completeness.",
        "revision/10q26/finngen_susie_*.tsv; ld_finngen_sisu42_pairs.json; coloc_bf_bf_AMDsusie_vs_stroke.tsv; ld_*.tsv; effective_N_region.tsv; tissue_eqtl/susie_10q26_*.tsv")

# ======================= S13 tissue eQTL + HTRA1 pQTL
B.sheet(13, "Multi-tissue eQTL colocalisation, top-eQTL SMR and HTRA1 pQTL MR for shared loci",
        "coloc.abf and coloc-SuSiE between tissue eQTL (GTEx v8, ROSMAP, BrainSeq, eQTLGen, EyeGEx retina) and AMD/stroke/WMH for seven shared regions; Wald-ratio SMR with the top cis-eQTL; key-variant eQTL effects; eQTL coverage; HTRA1 cis-pQTL MR (deCODE).",
        [("coloc.abf, tissue eQTL vs GWAS", rd("tissue_eqtl/coloc_tissue_eqtl_gwas.tsv")), ("coloc-SuSiE at 10q26", rd("tissue_eqtl/coloc_susie_10q26.tsv")),
         ("SMR using the top cis-eQTL (Wald ratio)", rd("tissue_eqtl/smr_tissue_topeqtl.tsv")),
         ("Key variants: tissue eQTL effects", rd("tissue_eqtl/key_variant_eqtl_lookup.tsv")),
         ("Retina/RPE/choroid single-nucleus eQTL look-up, 10q26 and 3p14 (Voigt et al.; FDR < 0.05 results only)", rd("revision3/eye_eqtl/sc_eqtl_10q26_3p14.tsv")),
         ("RPE HTRA1 eQTL: bounded conditional analysis of the two 10q26 haplotypes", rd("revision3/eye_eqtl/rpe_htra1_bound_conditional.tsv")),
         ("10q26 eQTLs conditional in both directions (retina, blood, arteries)", rd("revision5/eqtl_conditional_both_directions.tsv")),
         ("Retina and RPE-choroid bulk eQTL look-up (Orozco et al.; FDR <= 0.05 results only)", rd("revision3/eye_eqtl/bulk_eqtl_10q26.tsv")),
         ("eQTL coverage per region, gene and dataset", rd("tissue_eqtl/eqtl_dataset_gene_summary.tsv")),
         ("HTRA1 plasma pQTL MR (deCODE; instrument rs61871680)", rr("10q26/htra1_pqtl_mr_decode.tsv"))],
        "dir_expr_to_risk = sign of the effect of higher expression on disease risk; coloc_class = strong (PP.H4 > 0.8), moderate (0.5-0.8), distinct variants (PP.H3 > 0.8) or no evidence; b_allele = eQTL effect of the listed allele. The HTRA1 pQTL instrument is in low LD with both 10q26 signals (r2 <= 0.05).",
        "tissue_eqtl/; revision/10q26/htra1_pqtl_mr_decode.tsv")

# ======================= S14 3p14 PheWAS / coloc / HyPrColoc
B.sheet(14, "3p14 (ADAMTS9) haplotype: phenome-wide look-up, pairwise and multi-trait colocalisation, FinnGen credible sets and tissue scan",
        "PheWAS of rs4611812 and its haplotype (GWAS Catalog, FinnGen R12); pairwise coloc.abf with p12 sensitivity; HyPrColoc across AMD, aortic ADAMTS9 eQTL, stroke, cIMT and cardiometabolic traits; FinnGen in-sample coloc.bf_bf and credible-set overlap; Pan-UK Biobank replication with in-sample LD (AMD association not replicated); ADAMTS9 eQTL scan across 50 datasets; traits used.",
        [("PheWAS of the rs4611812 haplotype (curated; effects for the C allele)", rr("adamts9/phewas_rs4611812_haplotype.tsv")),
         ("GWAS Catalog associations in the region", rr("adamts9/gwascatalog_phewas.tsv")),
         ("FinnGen R12 PheWAS of rs4611812", rr("adamts9/finngen_phewas_rs4611812.tsv")),
         ("Pairwise coloc.abf (p12 = 1e-5, 5e-6, 1e-6)", rr("adamts9/coloc_pairwise_3p14.tsv")),
         ("HyPrColoc", rr("adamts9/hyprcoloc_3p14.tsv")),
         ("FinnGen in-sample coloc.bf_bf", rr("adamts9/coloc_bf_bf_finngen_insample.tsv")),
         ("Credible-set overlap", rr("adamts9/credible_set_overlap.tsv")), ("Credible-set members", rr("adamts9/credible_set_members.tsv")),
         ("Pan-UK Biobank replication: AMD effect of rs4611812-C and comparison with FinnGen/IAMDGC", rd("revision4/3p14_amd_replication.tsv")),
         ("Pan-UK Biobank phenotypes at 3p14", rd("revision4/3p14_panukb_phenotypes.tsv")),
         ("Pan-UK Biobank coloc.abf (in-sample LD)", rd("revision4/3p14_panukb_coloc.tsv")),
         ("Pan-UK Biobank SuSiE (in-sample LD)", rd("revision4/3p14_panukb_susie_summary.tsv")),
         ("ADAMTS9 eQTL tissue scan (50 testable datasets)", rr("adamts9/tissue_scan_ADAMTS9.tsv")),
         ("Traits used in locus analyses", rr("adamts9/traits_used.tsv"))],
        "beta_C_allele = effect of rs4611812-C (AMD-protective, aortic ADAMTS9-raising); posterior_prob = HyPrColoc posterior for the cluster; candidate_snp = most likely shared variant; analysis = trait set and prior (default, prior.c 0.05, reg/align 0.7); rank / p_bonf_datasets = rank of the dataset by minimum P and Bonferroni correction across 50 datasets.",
        "revision/adamts9/")

# ======================= S15 3p14 MVMR / heterogeneity / replication
def bp(f, lab):
    x = rr(f); x.insert(0, "blood_pressure_trait", lab); return x
B.sheet(15, "3p14 (ADAMTS9): mediation by cardiometabolic traits, subtype heterogeneity, coronary disease, replication and dry versus wet AMD",
        "Univariable and multivariable MR of T2D, WHRadjBMI and blood pressure on stroke (3p14 +/- 1 Mb excluded); expected versus observed stroke effects of rs4611812-C; heterogeneity between stroke subtypes (10q26 and 3p14 variants); stroke versus CAD; replication in MEGASTROKE, GIGASTROKE minus MEGASTROKE, FinnGen and East Asian GIGASTROKE; dry versus wet AMD.",
        [("Univariable MR, cardiometabolic -> stroke (SBP model)", bp("adamts9/mr_cardiometabolic_to_stroke.tsv", "SBP")),
         ("Univariable MR, cardiometabolic -> stroke (DBP model)", bp("adamts9/mr_cardiometabolic_to_stroke_DBP.tsv", "DBP")),
         ("MVMR, T2D + WHRadjBMI + SBP -> stroke", bp("adamts9/mvmr_cardiometabolic_to_stroke.tsv", "SBP")),
         ("MVMR, T2D + WHRadjBMI + DBP -> stroke", bp("adamts9/mvmr_cardiometabolic_to_stroke_DBP.tsv", "DBP")),
         ("Expected vs observed effect of rs4611812-C (SBP model)", rr("adamts9/expected_vs_observed_rs4611812.tsv")),
         ("Expected vs observed effect of rs4611812-C (DBP model)", rr("adamts9/expected_vs_observed_rs4611812_DBP.tsv")),
         ("Stroke subtype heterogeneity (10q26 and 3p14 variants)", rr("adamts9/subtype_heterogeneity.tsv")),
         ("Ischaemic stroke vs coronary disease", rr("adamts9/stroke_vs_cad_key_variants.tsv")),
         ("Replication of rs4611812-C", rr("adamts9/replication_rs4611812C.tsv")),
         ("Dry vs wet AMD", rr("adamts9/dry_vs_wet_amd.tsv"))],
        "theta = causal estimate per unit (T2D log-odds, WHRadjBMI SD, BP mmHg); exp_mvmr = stroke effect of rs4611812-C expected through the modelled risk factors; prop_explained_mvmr = exp_mvmr / obs_beta; residual_beta / residual_p = observed minus expected; comparison = subtype pair or three-subtype Q (shared-control correlation included); GIGASTROKE minus MEGASTROKE = approximate inverse-variance subtraction; rho (dry vs wet) = correlation of estimates from overlapping cases and shared controls.",
        "revision/adamts9/")

# ======================= S16 FDR tiers
fc = rr("mr/fdr_cited_results.tsv")
extra_layers = fa[fa["layer"].str.startswith("MR replication") | fa["layer"].str.startswith("MR to CSVD")].copy()
lay = fa.groupby(["layer", "tier"], as_index=False).agg(n_tests=("fdr", "size"), n_fdr_lt_0_05=("fdr", lambda s: int((s < 0.05).sum())))
B.sheet(16, "False discovery rate tiers for all results quoted in the manuscript",
        "One Benjamini-Hochberg family per analytical layer; every result quoted in the text with its layer, tier (primary / replication / exploratory), number of tests, FDR and status; FinnGen-free MEGASTROKE replication and CSVD-extension MR layers in full; number of tests and FDR-significant results per layer.",
        [("Results quoted in the manuscript", fc),
         ("MEGASTROKE replication and CSVD-extension MR layers (all tests)", extra_layers),
         ("Summary of FDR families", lay)],
        "tier = primary (defined after the initial analyses; not preregistered) / replication / exploratory; n_tests = size of the BH family; status = FDR < 0.05, nominal only (exploratory) or not significant. Full layer-level results (171,277 tests): results/revision/mr/fdr_all_layers.tsv.",
        "revision/mr/fdr_cited_results.tsv; revision/mr/fdr_all_layers.tsv")

# ======================= S17 SMR-HEIDI + pQTL
def nodrug(df):
    return df[[c for c in df.columns if 'drug' not in c.lower() and 'finan' not in c.lower()]]
sa = rd("smr_heidi/smr_heidi_all.tsv", low_memory=False)
sa_sig = sa[sa["p_fdr"] < 0.05].copy()
pa = rd("pqtl_mr/pqtl_mr_all.tsv"); pa_sig = pa[pa["fdr"] < 0.05].copy()
xq = (R / "pqtl_mr/xqtl_intersections.txt").read_text(encoding="utf-8", errors="ignore").strip().splitlines()
B.sheet(17, "Transcriptome-wide SMR-HEIDI (eQTLGen) and plasma protein cis-pQTL Mendelian randomization",
        "SMR counts per outcome; all probe-outcome pairs with FDR < 0.05; genes shared by AMD and stroke/WMH (including the 13 genes shared before HEIDI) with HEIDI and colocalisation; eQTL-GWAS colocalisation; pQTL MR results with FDR < 0.05, proteins shared by AMD and stroke/WMH, deCODE replication and xQTL intersections.",
        [("SMR counts per outcome", rd("smr_heidi/smr_heidi_counts.tsv")), ("SMR results with FDR < 0.05", nodrug(sa_sig)),
         ("Genes shared by AMD and stroke/WMH (FDR + HEIDI; relaxed criterion)", nodrug(rd("smr_heidi/druggable_shared_summary.tsv"))),
         ("Single-SNP SMR for the 82 shared candidate genes (before HEIDI)", rd("smr/smr_shared_genes.tsv")),
         ("eQTL-GWAS colocalisation for SMR candidates", rd("smr_heidi/coloc_eqtl_gwas.tsv")),
         ("pQTL MR results with FDR < 0.05 (UKB-PPP and deCODE)", pa_sig), ("Proteins shared by AMD and stroke/WMH", rd("pqtl_mr/pqtl_mr_shared.tsv")),
         ("Replication in deCODE", rd("pqtl_mr/pqtl_replication_decode.tsv")),
         ("Complement proteins: summary of pQTL MR across platforms with instrument annotation", rd("story/complement_pqtl_summary.tsv")), ("xQTL intersections", pd.DataFrame({"intersection": xq}))],
        "SMR robust = FDR < 0.05 within outcome and P_HEIDI >= 0.05; pQTL MR = Wald ratio with one cis variant per protein, OR per SD; replicated_decode = same direction and P < 0.05 in deCODE; region flags pleiotropic LD regions.",
        f"smr_heidi/ (full SMR table {len(sa):,} rows: smr_heidi_all.tsv); smr/smr_shared_genes.tsv; pqtl_mr/ (full pQTL table {len(pa):,} rows: pqtl_mr_all.tsv)")

# ======================= S18 single cell
sc_blocks = []
for ds in ("rpe_choroid", "retina_rpe_cowan", "brain_vascular_siletti", "brain_vascular_multiomics"):
    sc_blocks.append((f"ADAMTS9 and ADAMTS9-AS2 by cell type: {ds}", rd(f"adamts9/scrna/{ds}_celltype_expr.tsv")))
genes = ["HTRA1", "ARMS2", "ADAMTS9", "SH2B3", "ABO", "C4A", "C4B", "CFH", "C3"]
long = []
for ds in ("rpe_choroid", "retina_rpe_cowan", "brain_vascular_siletti", "brain_vascular_multiomics"):
    mn = pd.read_csv(R / f"scrna/{ds}/celltype_mean_expr.tsv", sep="\t").set_index("ct")
    fr = pd.read_csv(R / f"scrna/{ds}/celltype_fraction_expr.tsv", sep="\t").set_index("ct")
    for g in genes:
        if g in mn.columns:
            for ct in mn.index:
                long.append(dict(dataset=ds, gene=g, cell_type=ct, mean_log_expr=mn.at[ct, g], frac_expressing=fr.at[ct, g] if g in fr.columns else None))
long = pd.DataFrame(long)
B.sheet(18, "Single-cell expression of prioritised genes by cell type",
        "ADAMTS9 and ADAMTS9-AS2 mean expression and detection fraction per cell type in four CELLxGENE datasets (normal donors; author-provided annotations), dataset summary, and expression of other prioritised genes (HTRA1, ARMS2, SH2B3, ABO, C4A, C4B, CFH, C3) by cell type.",
        sc_blocks + [("ADAMTS9 / ADAMTS9-AS2 summary per dataset", rd("adamts9/scrna/summary.tsv")),
                     ("Prioritised genes by cell type (long format)", long)],
        "mean = mean log-normalised expression (log1p CP10k) in the cell type; frac = fraction of cells with detectable expression; n_cells = cells in the cell type; top3 = three cell types with the highest mean (mean/fraction).",
        "adamts9/scrna/*.tsv; scrna/*/celltype_mean_expr.tsv; scrna/*/celltype_fraction_expr.tsv")

ABBR = [("AMD", "age-related macular degeneration"), ("AS / IS / CES / LAS / SVS", "any / ischaemic / cardioembolic / large artery / small vessel stroke"),
        ("MEGA_*", "MEGASTROKE European outcomes (no FinnGen)"), ("EAS", "East Asian"), ("ICH", "intracerebral haemorrhage"),
        ("CSVD", "cerebral small vessel disease"), ("WMH", "white matter hyperintensities"), ("WM-MD / WM-FA", "white matter mean diffusivity / fractional anisotropy"),
        ("EPVS", "enlarged perivascular spaces (WM white matter, BG basal ganglia, HIP hippocampus)"), ("cIMT", "carotid intima-media thickness"),
        ("CAD / CHD", "coronary artery disease / coronary heart disease"), ("T2D", "type 2 diabetes"), ("WHRadjBMI", "waist-to-hip ratio adjusted for body mass index"),
        ("SBP / DBP", "systolic / diastolic blood pressure"), ("FG", "FinnGen"), ("IAMDGC", "International AMD Genomics Consortium"),
        ("A_dist / B_52 / C_orig", "instrument schemes: distance-defined loci (primary) / 52 conditionally independent IAMDGC variants / LD clumping"),
        ("MR", "Mendelian randomization"), ("IVW", "inverse-variance weighted"), ("MR-RAPS", "robust adjusted profile score"), ("MVMR", "multivariable MR"),
        ("TOST", "two one-sided tests of equivalence"), ("MDE", "minimum detectable effect"), ("LDSC", "linkage disequilibrium score regression"),
        ("LAVA", "Local Analysis of [co]Variant Association"), ("SMR / HEIDI", "summary-data-based MR / heterogeneity in dependent instruments"),
        ("PP.H3 / PP.H4", "posterior probability of distinct / shared causal variants"), ("PIP", "posterior inclusion probability"),
        ("COJO", "conditional and joint analysis (GCTA)"), ("SuSiE", "sum of single effects fine-mapping"), ("eQTL / pQTL", "expression / protein quantitative trait locus"),
        ("DLPFC", "dorsolateral prefrontal cortex"), ("RPE", "retinal pigment epithelium"), ("OR", "odds ratio"), ("CI", "confidence interval"),
        ("FDR", "false discovery rate (Benjamini-Hochberg)"), ("Neff", "effective sample size")]
B.save(ABBR)
(Path(f"{PROJ_ROOT}/submission/revision/build") / "Additional_file_2_build_log.txt").write_text("\n".join(LOG), encoding="utf-8")
print("\n".join(LOG)); print("saved", OUT)
