"""
Sets up CPHI-DRAGEN-anno analysis directories and downloads HPO terms from Phenotips.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import glob
import logging
import pandas as pd 
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


logger = logging.getLogger("DRAGEN_setup")

REPO_ROOT_DEFAULT_CPHI_DRAGEN = Path.home() / "CPHI-DRAGEN-anno"
DEFAULT_CREDS = "PT_credentials.csv"

BASE = Path("/hpf/largeprojects/tgnode/sandbox/mcouse_analysis")
HPO_DIR = BASE / "HPO"
PED_DIR = BASE / "pedigrees"
ANALYSES_BASE = BASE / "analyses" 
PCHSEQ_DIR = Path("/hpf/projects/PCHSeq/")
PROJECT_DICT = {"sickkidsseq": "SickKidsSeq", "genoderm": "SkinGene"}


@dataclass(frozen=True)
class AnalysisRow:
    family: str
    sequence_id: str
    project_id: str
    sample_type: str
    family_pchseq: str
    lims: str

def _configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

def _run(cmd: list[str], *, cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    logger.debug("Running command: %s (cwd=%s)", " ".join(cmd), cwd)
    try:
        result = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            check=True,
            text=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as e:
        logger.error("Command failed (exit %s): %s", e.returncode, " ".join(cmd))
        if e.stdout:
            logger.error("stdout:\n%s", e.stdout)
        if e.stderr:
            logger.error("stderr:\n%s", e.stderr)
        raise
    if result.stdout:
        logger.debug("stdout:\n%s", result.stdout)
    if result.stderr:
        logger.debug("stderr:\n%s", result.stderr)
    return result

def _strip_cr(s: str) -> str:
    return s.replace("\r", "").strip()

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def parse_analysis_tsv(path: Path) -> list[AnalysisRow]:
    logger.info("Parsing analysis TSV: %s", path)
    rows: list[AnalysisRow] = []
    analysis_df = pd.read_csv(path, sep="\t")
    for index, row in analysis_df.iterrows():
        family = row["Family_ID"]
        family_pchseq = row["Family_ID_PCHseq"]
        sequence_id = row["Sequence_ID"]
        project_id = row["TG_ID"]
        sample_type = row["Sample_type"]
        lims = row["LIMS"]
        rows.append(AnalysisRow(family=family, family_pchseq=family_pchseq, sequence_id=sequence_id, project_id=project_id, sample_type=sample_type, lims=lims))

    logger.info("Parsed %d row(s) from %s", len(rows), path.name)
    return rows


def normalize_family_id(family: str) -> str:
    """
    Normalize family IDs so directory names are consistent across runs.
    - Drop anything after a '.' (some inputs include suffixes)
    - Remove '_' and '-' characters
    """
    return family.split(".", 1)[0].replace("_", "").replace("-", "")

def init_existing_family_dirs(analysis_rows: list[AnalysisRow], analysis_dir: Path, today: str) -> None:
    """
    Set up analysis directories for families that already exist.
    """
    logger.info("Resetting samples.tsv for any pre-existing family dirs under %s", analysis_dir)
    seen: set[str] = set()
    reset_count = 0
    for r in analysis_rows:
        family_norm = normalize_family_id(r.family)
        if not family_norm or family_norm in seen:
            continue
        seen.add(family_norm)
        family_dir = analysis_dir / family_norm / f"DRAGEN_{today}"
        if family_dir.is_dir():
            samples_tsv = family_dir / "samples.tsv"
            samples_tsv.write_text("sample\tDRAGEN_results_dir\n")
            logger.debug("Reset samples.tsv in existing dir: %s", family_dir)
            reset_count += 1
    logger.info("Reset samples.tsv for %d existing family dir(s)", reset_count)

def rewrite_config_yaml(config_path: Path, cphi: bool, *, family: str, hpo: Optional[Path], ped: Optional[Path]) -> None:
    logger.debug("Rewriting config %s (family=%s, hpo=%s, ped=%s)", config_path, family, hpo, ped)
    txt = config_path.read_text()
    txt = txt.replace("FAM-000000", family)
    txt = txt.replace("~/", "/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/tools/")
    if cphi:
        txt = txt.replace("cphi: false", "cphi: true")
    else:
        txt = txt.replace(' dragen_output_schema: ""', 'dragen_output_schema: "modified"')
    if hpo is not None:
        txt = txt.replace('hpo: ""', f'hpo: "{hpo}"')
    else:
        logger.warning("No HPO file substituted in %s", config_path)
    if ped is not None:
        txt = txt.replace('ped: ""', f'ped: "{ped}"')
    else:
        logger.warning("No pedigree file substituted in %s", config_path)
    config_path.write_text(txt)

def find_hpo(project: str, family: str, family_norm: str) -> Optional[Path]:
    base = HPO_DIR / project
    if ("DSK" in family_norm) or ("GYM" in family_norm) or ("SKS" in family_norm):
        matches = sorted(base.glob(f"{family_norm}*"), key=lambda p: p.stat().st_mtime, reverse=True)
        hit = matches[0] if matches else None
    elif "GD" in family:
        matches = sorted(base.glob(f"{family}*"), key=lambda p: p.stat().st_mtime, reverse=True)
        hit = matches[0] if matches else None
    else:
        hit = None
    if hit is None:
        logger.warning("No HPO file found for family=%s (normalized=%s) under %s", family, family_norm, base)
    else:
        logger.debug("HPO file for family=%s: %s", family, hit)
    return hit

def find_pedigree(DRAGEN_joint_geno_dir: Path, family_pchseq: str) -> Optional[Path]:
    ped = DRAGEN_joint_geno_dir / f"{family_pchseq}.ped"
    if not ped.exists():
        logger.warning("Expected pedigree not found: %s", ped)
    return ped

def find_pedigree_nonCPHI(project: str, family_norm: str, family: str) -> Optional[Path]:
    base = PED_DIR / project
    if ("DSK" in family) or ("GYM" in family) or ("SKS" in family):
        matches = sorted(base.glob(f"{family}*"), key=lambda p: p.stat().st_mtime, reverse=True)
        hit = matches[0] if matches else None
    elif "GD" in family:
        matches = sorted(base.glob(f"{family}*"), key=lambda p: p.stat().st_mtime, reverse=True)
        hit = matches[0] if matches else None
    else:
        hit = None
    if hit is None:
        logger.warning("No pedigree file found for family=%s (normalized=%s) under %s", family, family_norm, base)
    else:
        logger.debug("Pedigree file for family=%s: %s", family, hit)
    return hit

def setup_family_once(
    *,
    family: str,
    project: str,
    family_norm: str,
    family_pchseq: str,
    family_dir: Path,
    lims: str,
    cphi: bool,
    cphi_dragen_anno: Path,
    today: str,
) -> None:
    """
    Ensure family dir exists and has config + units.tsv + samples.tsv, plus cnv/str dirs.
    Returns (sequence_variant_vcf, SV_vcf, CNV_vcf, repeat_VCF_dir, VNTR_VCF_dir).
    """
    if cphi:
        DRAGEN_joint_geno_dir = PCHSEQ_DIR / PROJECT_DICT[project] / f"{lims}_family" / family_pchseq / "output" 
        sequence_variant_vcf = DRAGEN_joint_geno_dir / f"{family_pchseq}.hard-filtered.vcf.gz"
        SV_vcf = DRAGEN_joint_geno_dir / f"{family_pchseq}.sv.vcf.gz"
        CNV_vcf = DRAGEN_joint_geno_dir / f"{family_pchseq}.cnv.vcf.gz"
    else:
        sequence_variant_vcf = glob.glob(f"/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/{project}/{lims}/FAM*{family}*hard-filtered.vcf.gz")[0]
        SV_vcf = glob.glob(f"/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/{project}/{lims}/FAM*{family}*sv.with-inv.vcf.gz")[0]
        CNV_vcf = glob.glob(f"/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/{project}/{lims}/FAM*{family}*cnv.vcf.gz")[0]

    if family_dir.exists():
        logger.info("Family dir already exists, skipping initial setup: %s", family_dir)
        return

    logger.info("Setting up family dir: %s (family=%s, family_pchseq=%s)", family_dir, family, family_pchseq)
    ensure_dir(family_dir)

    for src in [
        cphi_dragen_anno / "config" / "config.yaml",
        cphi_dragen_anno / "workflow" / "CPHI_DRAGEN_anno.sh",
    ]:
        if not src.exists():
            logger.error("Missing pipeline file: %s", src)
            raise FileNotFoundError(f"Missing pipeline file: {src}")
        logger.debug("Copying %s -> %s", src, family_dir / src.name)
        shutil.copy2(src, family_dir / src.name)

    config_path = family_dir / "config.yaml"
    hpo = find_hpo(project, family, family_norm)
    if cphi:
        ped = find_pedigree(DRAGEN_joint_geno_dir, family_pchseq)
    else:
        ped = find_pedigree_nonCPHI(project, family_norm, family)
    logger.debug("Copying pedigree %s -> %s", ped, family_dir / f"{family_pchseq}.ped")
    shutil.copy2(ped, family_dir / f"{family_pchseq}.ped")
    if cphi:
        rewrite_config_yaml(config_path, cphi, family=family_pchseq, hpo=hpo, ped=ped)
    else:
        rewrite_config_yaml(config_path, cphi, family=family_norm, hpo=hpo, ped=ped)

    (family_dir / "samples.tsv").write_text("sample\tDRAGEN_results_dir\n")

    units_tsv = family_dir / "units.tsv"
    units_tsv.write_text("family\tsmall_variant_vcf\tSV_vcf\tCNV_vcf\n")
    if cphi:
        with units_tsv.open("a") as out:
            out.write(f"{family_pchseq}\t{sequence_variant_vcf}\t{SV_vcf}\t{CNV_vcf}\n")
        logger.debug("Wrote units.tsv for %s", family_pchseq)

    else:
        with units_tsv.open("a") as out:
            out.write(f"{family_norm}\t{sequence_variant_vcf}\t{SV_vcf}\t{CNV_vcf}\n")
        logger.debug("Wrote units.tsv for %s", family_norm)


def submit_cphi_dragen_anno_slurm(family_dir: Path) -> None:
    """
    Submit the CPHI_DRAGEN_anno Slurm script from the family analysis directory.
    """
    script = family_dir / "CPHI_DRAGEN_anno.sh"
    if not script.is_file():
        logger.error("Expected Slurm script not found: %s", script)
        raise FileNotFoundError(f"Expected Slurm script not found: {script}")
    logger.info("Submitting Slurm job: sbatch %s (cwd=%s)", script, family_dir)
    _run(["sbatch", str(script)], cwd=family_dir)


def add_sample_inputs(
    *,
    family_dir: Path,
    family_pchseq: str,
    project: str,
    sequence_id: str,
    lims: str,
    cphi: bool,
) -> None:
    if cphi:
        dragen_results_dir_sample = PCHSEQ_DIR / PROJECT_DICT[project] / lims / f"{sequence_id}" 
    else:
        dragen_results_dir_sample = f"/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/{project}/{lims}/"
    if not Path(dragen_results_dir_sample).exists():
        raise FileNotFoundError(f"DRAGEN results dir does not exist: {dragen_results_dir_sample}")
    logger.info("Adding sample %s -> %s/samples.tsv", sequence_id, family_dir)
    with (family_dir / "samples.tsv").open("a") as out:
        out.write(f"{sequence_id}\t{dragen_results_dir_sample}\n")



def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="Set up CPHI-DRAGEN-anno analysis directories for PCHseq DRAGEN inputs.")
    ap.add_argument("--analyses", type=Path, help="Path to sample metadata TSV")
    ap.add_argument("--project", help="Project ID, e.g. DECODER")
    ap.add_argument("--cphi", default="True", help="True if CPHI/PCHseq sample, otherwise False")
    ap.add_argument("--creds", default=DEFAULT_CREDS, help="Phenotips credentials CSV (default: PT_credentials.csv)")
    ap.add_argument("--cphi-dragen-anno", dest="cphi_dragen_anno", type=Path, default=REPO_ROOT_DEFAULT_CPHI_DRAGEN, help="Path to CPHI-DRAGEN-anno repo (default: ~/CPHI-DRAGEN-anno)")
    ap.add_argument("--today", default=None, help="Override date stamp (YYYY-MM-DD). Default: today.")
    ap.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging verbosity (default: INFO).",
    )
    args = ap.parse_args(argv)

    _configure_logging(args.log_level)

    logger.info("Starting DRAGEN_setup (project=%s, analyses=%s)", args.project, args.analyses)

    if not args.analyses.is_file():
        logger.error("Sample file does not exist: %s", args.analyses)
        raise SystemExit(f"Error: sample file does not exist: {args.analyses}")

    today = args.today or _dt.date.today().isoformat()
    logger.info("Using date stamp: %s", today)

    analysis_dir = ANALYSES_BASE / args.project
    ensure_dir(analysis_dir)
    logger.info("Analysis base dir: %s", analysis_dir)

    cphi = args.cphi.lower() == "true"
    logger.info("Using CPHI: %s", cphi)

    rows = parse_analysis_tsv(args.analyses)

    logger.info("Downloading HPO terms and pedigrees from Phenotips")
    _run(
        [
            "python3",
            "get_HPO_pedigree_genome_clinic.py",
            "-sample_sheet",
            str(args.analyses),
            "-credentials",
            args.creds,
            "-project",
            args.project,
            "-rename",
            "False",
        ]
    )

    init_existing_family_dirs(rows, analysis_dir, today)

    last_row_index_by_family: dict[str, int] = {}
    for idx, r in enumerate(rows):
        last_row_index_by_family[normalize_family_id(r.family)] = idx

    logger.info("Processing %d analysis row(s)", len(rows))
    for idx, r in enumerate(rows):
        family = r.family
        sequence_id = _strip_cr(r.sequence_id)

        family_norm = normalize_family_id(family)
        family_pchseq = r.family_pchseq
        family_dir = analysis_dir / family_norm / f"DRAGEN_{today}"

        logger.info("Processing family=%s (norm=%s, pchseq=%s, sample=%s)", family, family_norm, family_pchseq, sequence_id)

        try:
            setup_family_once(
                family=family,
                lims=r.lims,
                project=args.project,
                family_norm=family_norm,
                family_pchseq=family_pchseq,
                family_dir=family_dir,
                cphi=cphi,
                cphi_dragen_anno=args.cphi_dragen_anno,
                today=today,
            )

            add_sample_inputs(
                lims=r.lims,
                family_dir=family_dir,
                family_pchseq=family_pchseq,
                project=args.project,
                sequence_id=sequence_id,
                cphi=cphi,
            )

            if idx == last_row_index_by_family[family_norm]:
                submit_cphi_dragen_anno_slurm(family_dir)
        except Exception:
            logger.exception("Failed to set up family=%s (sample=%s)", family, sequence_id)
            raise

    logger.info("DRAGEN_setup complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

