"""
Sets up CPHI-DRAGEN-anno analysis directories and downloads HPO terms from Phenotips.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import glob
import pandas as pd 
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


REPO_ROOT_DEFAULT_CPHI_DRAGEN = Path.home() / "CPHI-DRAGEN-anno"
DEFAULT_CREDS = "PT_credentials.csv"

BASE = Path("/hpf/largeprojects/tgnode/sandbox/mcouse_analysis")
HPO_DIR = BASE / "HPO"
PED_DIR = BASE / "pedigrees"
FILES_FROM_IRODS = BASE / "files_from_irods"
ANALYSES_BASE = BASE / "analyses" 

@dataclass(frozen=True)
class AnalysisRow:
    family: str
    sequence_id: str
    project_id: str
    sample_type: str
    family_pchseq: str

def _run(cmd: list[str], *, cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        check=True,
        text=True,
        capture_output=True,
    ) 

def _strip_cr(s: str) -> str:
    return s.replace("\r", "").strip()

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def parse_analysis_tsv(path: Path) -> list[AnalysisRow]:
    rows: list[AnalysisRow] = []
    analysis_df = pd.read_csv(path, sep="\t")
    for index, row in analysis_df.iterrows():
        family = row["Family_ID"]
        family_pchseq = row["Family_ID_PCHseq"]
        sequence_id = row["Sequence_ID"]
        project_id = row["TG_ID"]
        sample_type = row["Sample_type"]
        rows.append(AnalysisRow(family=family, family_pchseq=family_pchseq, sequence_id=sequence_id, project_id=project_id, sample_type=sample_type))

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
    seen: set[str] = set()
    for r in analysis_rows:
        family_norm = normalize_family_id(r.family)
        if not family_norm or family_norm in seen:
            continue
        seen.add(family_norm)
        family_dir = analysis_dir / family_norm / f"DRAGEN_{today}"
        if family_dir.is_dir():
            samples_tsv = family_dir / "samples.tsv"
            samples_tsv.write_text("sample\tDRAGEN_results_dir\n")

def rewrite_config_yaml(config_path: Path, *, family_pchseq: str, hpo: Optional[Path], ped: Optional[Path]) -> None:
    txt = config_path.read_text()
    txt = txt.replace("FAM-000000", family_pchseq)
    if hpo is not None:
        txt = txt.replace('hpo: ""', f'hpo: "{hpo}"')
    if ped is not None:
        txt = txt.replace('ped: ""', f'ped: "{ped}"')
    config_path.write_text(txt)

def find_hpo(project: str, family: str, family_norm: str) -> Optional[Path]:
    base = HPO_DIR / project
    if ("DSK" in family_norm) or ("GYM" in family_norm):
        patt = str(base / f"{family_norm}*")
        matches = sorted(base.glob(f"{family_norm}*"), key=lambda p: p.stat().st_mtime, reverse=True)
        return matches[0] if matches else None
    if "GD" in family:
        matches = sorted(base.glob(f"{family}*"), key=lambda p: p.stat().st_mtime, reverse=True)
        return matches[0] if matches else None
    return None

def find_pedigree(DRAGEN_joint_geno_dir: Path, family_pchseq: str) -> Optional[Path]:
    return DRAGEN_joint_geno_dir / f"{family_pchseq}.ped"

def setup_family_once(
    *,
    family: str,
    project: str,
    family_norm: str,
    family_pchseq: str,
    family_dir: Path,
    cphi_dragen_anno: Path,
    today: str,
) -> None:
    """
    Ensure family dir exists and has config + units.tsv + samples.tsv, plus cnv/str dirs.
    Returns (sequence_variant_vcf, SV_vcf, CNV_vcf, repeat_VCF_dir, VNTR_VCF_dir).
    """
    DRAGEN_joint_geno_dir = FILES_FROM_IRODS / project / "PCHseq" / family_pchseq / "output" 
    sequence_variant_vcf = DRAGEN_joint_geno_dir / f"{family_pchseq}.hard-filtered.vcf.gz"
    SV_vcf = DRAGEN_joint_geno_dir / f"{family_pchseq}.sv.vcf.gz"
    CNV_vcf = DRAGEN_joint_geno_dir / f"{family_pchseq}.cnv.vcf.gz"

    if not family_dir.exists():
        ensure_dir(family_dir)

        # Copy pipeline files
        for src in [
            cphi_dragen_anno / "config" / "config.yaml",
            cphi_dragen_anno / "workflow" / "CPHI_DRAGEN_anno.sh",
        ]:
            if not src.exists():
                raise FileNotFoundError(f"Missing pipeline file: {src}")
            shutil.copy2(src, family_dir / src.name)

        config_path = family_dir / "config.yaml"
        hpo = find_hpo(project, family, family_norm)
        ped = find_pedigree(DRAGEN_joint_geno_dir, family_pchseq)
        # Copy pedigree to family dir
        shutil.copy2(ped, family_dir / f"{family_pchseq}.ped")
        rewrite_config_yaml(config_path, family_pchseq=family_pchseq, hpo=hpo, ped=family_dir / f"{family_pchseq}.ped")

        # Create samples.tsv / units.tsv
        (family_dir / "samples.tsv").write_text("sample\tDRAGEN_results_dir\n")

        units_tsv = family_dir / "units.tsv"
        units_tsv.write_text("family\tsequence_variant_vcf\tSV_vcf\tCNV_vcf\n")
        with units_tsv.open("a") as out:
            out.write(f"{family_pchseq}\t{sequence_variant_vcf}\t{SV_vcf}\t{CNV_vcf}\n")


def add_sample_inputs(
    *,
    family_dir: Path,
    family_pchseq: str,
    project: str,
    sequence_id: str,
) -> None:
    # samples.tsv
    dragen_results_dir_sample = FILES_FROM_IRODS / project / "PCHseq" / sequence_id 
    ensure_dir(dragen_results_dir_sample)
    with (family_dir / "samples.tsv").open("a") as out:
        out.write(f"{sequence_id}\t{dragen_results_dir_sample}\n")



def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="Set up crg2-pacbio analysis directories for PacBio inputs.")
    ap.add_argument("--analyses", type=Path, help="Path to sample metadata TSV")
    ap.add_argument("--project", help="Project ID, e.g. DECODER")
    ap.add_argument("--creds", default=DEFAULT_CREDS, help="Phenotips credentials CSV (default: PT_credentials.csv)")
    ap.add_argument("--cphi-dragen-anno", dest="cphi_dragen_anno", type=Path, default=REPO_ROOT_DEFAULT_CPHI_DRAGEN, help="Path to CPHI-DRAGEN-anno repo (default: ~/CPHI-DRAGEN-anno)")
    ap.add_argument("--today", default=None, help="Override date stamp (YYYY-MM-DD). Default: today.")
    args = ap.parse_args(argv)

    if not args.analyses.is_file():
        raise SystemExit(f"Error: sample file does not exist: {args.analyses}")

    today = args.today or _dt.date.today().isoformat()

    analysis_dir = ANALYSES_BASE / args.project
    ensure_dir(analysis_dir)

    rows = parse_analysis_tsv(args.analyses)

    # Download HPO + pedigrees from Phenotips (we will use pedigree from DRAGEN results directory instead)
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
        ]
    )

    # If analysis dir already exists for family, reset samples.tsv
    init_existing_family_dirs(rows, analysis_dir, today)

    for r in rows:
        family = r.family
        sequence_id = _strip_cr(r.sequence_id)

        family_norm = normalize_family_id(family)
        family_pchseq = r.family_pchseq
        family_dir = analysis_dir / family_norm / f"DRAGEN_{today}"

        setup_family_once(
            family=family,
            project=args.project,
            family_norm=family_norm,
            family_pchseq=family_pchseq,
            family_dir=family_dir,
            cphi_dragen_anno=args.cphi_dragen_anno,
            today=today,
        )

        add_sample_inputs(
            family_dir=family_dir,
            family_pchseq=family_pchseq,
            project=args.project,
            sequence_id=sequence_id,
        )

    # validate_pedigrees(rows, args.analyses, args.project)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

