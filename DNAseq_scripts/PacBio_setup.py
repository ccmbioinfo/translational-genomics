#!/usr/bin/env python3
"""
Python rewrite of PacBio_setup.sh

Sets up crg2-pacbio analysis directories, downloads HPO terms and pedigrees from Phenotips,
copies per-sample inputs, rewrites sample IDs in VCFs using bcftools, and validates that
samples in the analysis TSV are present in the Phenotips pedigree.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import glob
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Union


REPO_ROOT_DEFAULT_CRG2_PACBIO = Path.home() / "crg2-pacbio"
DEFAULT_CREDS = "PT_credentials.csv"

BASE = Path("/hpf/largeprojects/tgnode/sandbox/mcouse_analysis")
HPO_DIR = BASE / "HPO"
PED_DIR = BASE / "pedigrees"
FILES_FROM_IRODS = BASE / "files_from_irods"
ANALYSES_BASE = BASE / "analyses/test"


@dataclass(frozen=True)
class AnalysisRow:
    family: str
    sequence_id: str
    project_id_raw: str
    sample_type: str

    @property
    def family_is_header(self) -> bool:
        return self.family == "Family_ID"


def _run(cmd: list[str], *, cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        check=True,
        text=True,
        capture_output=True,
    )


def _which_or_die(exe: str) -> None:
    if shutil.which(exe) is None:
        raise SystemExit(f"Error: required executable not found on PATH: {exe}")


def _latest_matching(glob_pat: str) -> Optional[Path]:
    matches = [Path(p) for p in glob.glob(glob_pat)]
    if not matches:
        return None
    matches.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return matches[0]


_PathLikeOrGlob = Union[Path, str]


def _is_glob_pattern(s: str) -> bool:
    # Minimal glob metachar check (* ? [..])
    return any(ch in s for ch in ("*", "?", "["))


def _glob_latest(paths_or_globs: Iterable[_PathLikeOrGlob]) -> Optional[Path]:
    """
    Return the most-recently-modified match among:
      - literal paths that exist
      - glob patterns (e.g. "/path/to/*.vcf.gz")

    This fixes the bash->python mismatch where wildcard strings were previously treated
    as literal paths and dropped by `.exists()`.
    """
    expanded: list[Path] = []
    for item in paths_or_globs:
        s = str(item)
        if _is_glob_pattern(s):
            expanded.extend(Path(p) for p in glob.glob(s))
        else:
            p = Path(s)
            if p.exists():
                expanded.append(p)

    expanded = [p for p in expanded if p.exists()]
    if not expanded:
        return None
    expanded.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return expanded[0]


def _strip_cr(s: str) -> str:
    return s.replace("\r", "").strip()


def normalize_project_id(family: str, project_id: str) -> str:
    """
    Mirrors the bash normalization:
      - DSK/GYM: tr -d '_' | tr '.' '_' | tr -d '\r'
      - GD: sed 's/GD-/GD/' | tr '-' '_' | tr -d '\r'
    """
    project_id = _strip_cr(project_id)
    if ("DSK" in family) or ("GYM" in family):
        return project_id.replace("_", "").replace(".", "_")
    if "GD" in family:
        if project_id.startswith("GD-"):
            project_id = "GD" + project_id[len("GD-") :]
        return project_id.replace("-", "_")
    raise ValueError(f"Family {family} not found / unsupported")


def project_family_from_project_id(project_id: str) -> str:
    return project_id.split("_", 1)[0]


def project_sample_from_project_id(project_id: str) -> str:
    parts = project_id.split("_", 2)
    if len(parts) < 2:
        raise ValueError(f"Unexpected project_id format (expected *_*): {project_id}")
    return parts[1].strip().replace("\r", "").replace(" ", "")


def tcag_family_if_needed(family: str) -> Optional[str]:
    if "DSK" in family:
        return family.replace("DSK_", "DSK_DNA_", 1)
    if "GYM" in family:
        return family.replace("GYM_", "GYM_DNA_", 1)
    return None


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def copy_with_sidecars(src_vcfgz: Path, dest_dir: Path) -> Path:
    """
    Copy src_vcfgz plus any sidecar files matching src_vcfgz* into dest_dir.
    Returns destination .vcf.gz path.
    """
    ensure_dir(dest_dir)
    copied_vcfgz: Optional[Path] = None
    for f in sorted(src_vcfgz.parent.glob(src_vcfgz.name + "*")):
        if f.is_file():
            dest = dest_dir / f.name
            shutil.copy2(f, dest)
            if dest.name.endswith(".vcf.gz"):
                copied_vcfgz = dest
    if copied_vcfgz is None:
        raise FileNotFoundError(f"Did not copy a .vcf.gz for {src_vcfgz}")
    return copied_vcfgz


def rewrite_config_yaml(config_path: Path, *, project_family: str, hpo: Optional[Path], ped: Optional[Path]) -> None:
    txt = config_path.read_text()
    txt = txt.replace("NA12878", project_family)
    if hpo is not None:
        txt = txt.replace('hpo: ""', f'hpo: "{hpo}"')
    if ped is not None:
        txt = txt.replace('ped: ""', f'ped: "{ped}"')
    txt = txt.replace("c4r: True", "c4r: False")
    config_path.write_text(txt)


def pick_deepvariant(project: str, family: str, sequence_id: str) -> Path:
    base = FILES_FROM_IRODS / project
    sequence_fam = sequence_id.split("_")[0] # may need this for DECODER samples sequenced under C4R IDs, e.g. 1741_SK0123
    candidates = [
        base / f"{family}-cohort.joint.GRCh38.small_variants.phased.vcf.gz",
        base / f"{family}.joint.GRCh38.small_variants.phased.vcf.gz", # TCAG does not always consistently name joint-genotyped VCFs
        base / f"{sequence_id}.GRCh38.small_variants.phased.vcf.gz", # singleton sample, so no joint-genotyped VCF
        base / f"{family}-fam.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz", # file naming for older PacBio pipeline runs
        base / f"{sequence_fam}-fam.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz", # file naming for older PacBio pipeline runs
    ]
    tcag = tcag_family_if_needed(family)
    if tcag:
        candidates.append(base / f"{tcag}-cohort.joint.GRCh38.small_variants.phased.vcf.gz")
    picked = _glob_latest(candidates)
    if picked is None:
        raise FileNotFoundError(f"Could not find deepvariant VCF for family={family} sequence_id={sequence_id} project={project}")
    return picked


def pick_sv(project: str, family: str, sequence_id: str) -> Path:
    base = FILES_FROM_IRODS / project
    sequence_fam = sequence_id.split("_")[0] # may need this for DECODER samples sequenced under C4R IDs, e.g. 1741_SK0123
    candidates = [
        base / f"{family}-cohort.joint.GRCh38.structural_variants.phased.vcf.gz",
        base / f"{family}.joint.GRCh38.structural_variants.phased.vcf.gz",
        base / f"{sequence_id}.GRCh38.structural_variants.phased.vcf.gz", # singleton sample, so no joint-genotyped VCF
        base / f"{family}-fam.joint.GRCh38.pbsv.phased.vcf.gz", # file naming for older PacBio pipeline runs
        base / f"{sequence_fam}-fam.joint.GRCh38.pbsv.phased.vcf.gz", # file naming for older PacBio pipeline runs
    ]
    tcag = tcag_family_if_needed(family)
    if tcag:
        candidates.append(base / f"{tcag}-cohort.joint.GRCh38.structural_variants.phased.vcf.gz")
    picked = _glob_latest(candidates)
    if picked is None:
        raise FileNotFoundError(f"Could not find structural variants VCF for family={family} sequence_id={sequence_id} project={project}")
    return picked

def find_hpo(project: str, family: str, project_family: str) -> Optional[Path]:
    base = HPO_DIR / project
    if ("DSK" in family) or ("GYM" in family):
        patt = str(base / f"{project_family}*")
        matches = sorted(base.glob(f"{project_family}*"), key=lambda p: p.stat().st_mtime, reverse=True)
        return matches[0] if matches else None
    if "GD" in family:
        matches = sorted(base.glob(f"{family}*"), key=lambda p: p.stat().st_mtime, reverse=True)
        return matches[0] if matches else None
    return None


def find_pedigree(project: str, family: str, project_family: str) -> Optional[Path]:
    base = PED_DIR / project
    matches = sorted(base.glob(f"{project_family}*"), key=lambda p: p.stat().st_mtime, reverse=True)
    if matches:
        return matches[0]
    matches = sorted(base.glob(f"{family}*"), key=lambda p: p.stat().st_mtime, reverse=True)
    return matches[0] if matches else None


def bcftools_reheader_inplace(vcfgz: Path, mapping_file: Path) -> None:
    """
    Equivalent to:
      bcftools reheader -s mapping vcfgz | bcftools view -Oz -o vcfgz.vcf.gz
      mv vcfgz.vcf.gz vcfgz
    but implemented as a single shell pipeline to avoid Python-level pipe/file-object issues.
    """
    out_tmp = vcfgz.with_suffix(vcfgz.suffix + ".vcf.gz")  # ".vcf.gz" appended
    cmd = (
        f"bcftools reheader -s {mapping_file} {vcfgz} "
        f"| bcftools view -Oz -o {out_tmp}"
    )
    subprocess.run(cmd, shell=True, check=True)
    out_tmp.replace(vcfgz)


def bcftools_query_sample(vcfgz: Path) -> str:
    cp = _run(["bcftools", "query", "-l", str(vcfgz)])
    line = (cp.stdout or "").splitlines()
    if not line:
        raise RuntimeError(f"bcftools query -l returned no samples for {vcfgz}")
    return line[0].strip()


def parse_analysis_tsv(path: Path) -> list[AnalysisRow]:
    rows: list[AnalysisRow] = []
    with path.open(newline="") as f:
        rdr = csv.reader(f, delimiter="\t")
        for parts in rdr:
            if not parts:
                continue
            # Bash reads 10 columns but only uses family, sequence_id, project_id, sample_type
            family = _strip_cr(parts[0]) if len(parts) > 0 else ""
            sequence_id = _strip_cr(parts[1]) if len(parts) > 1 else ""
            project_id = _strip_cr(parts[2]) if len(parts) > 2 else ""
            sample_type = _strip_cr(parts[3]) if len(parts) > 3 else ""
            rows.append(AnalysisRow(family=family, sequence_id=sequence_id, project_id_raw=project_id, sample_type=sample_type))
    return rows


def init_existing_family_dirs(analysis_rows: list[AnalysisRow], analysis_dir: Path, today: str) -> None:
    """
    Mirrors:
      for project_family in `awk '{print $3}' analyses | cut -d '.' -f1 | tr -d '_' | uniq`
        if dir exists: recreate samples.tsv with header
    """
    seen: set[str] = set()
    for r in analysis_rows:
        if r.family_is_header:
            continue
        project_family = r.project_id_raw.split(".", 1)[0].replace("_", "")
        if not project_family or project_family in seen:
            continue
        seen.add(project_family)
        family_dir = analysis_dir / project_family / f"PacBio_{today}"
        if family_dir.is_dir():
            samples_tsv = family_dir / "samples.tsv"
            samples_tsv.write_text("sample\tbam\tcase_or_control\n")


def setup_family_once(
    *,
    family: str,
    sequence_id: str,
    project: str,
    project_id_norm: str,
    project_family: str,
    family_dir: Path,
    crg2_pacbio: Path,
    today: str,
) -> tuple[Path, Path]:
    """
    Ensure family dir exists and has config + units.tsv + samples.tsv, plus cnv/trgt dirs.
    Returns (deepvariant_vcfgz, sv_vcfgz).
    """
    deepvariant = pick_deepvariant(project, family, sequence_id)
    sv = pick_sv(project, family, sequence_id)

    if not family_dir.exists():
        ensure_dir(family_dir)

        # Copy pipeline files
        for src in [
            crg2_pacbio / "config.yaml",
            crg2_pacbio / "crg2-pacbio.sh",
            crg2_pacbio / "slurm_profile" / "slurm-config.yaml",
        ]:
            if not src.exists():
                raise FileNotFoundError(f"Missing pipeline file: {src}")
            shutil.copy2(src, family_dir / src.name)

        config_path = family_dir / "config.yaml"
        hpo = find_hpo(project, family, project_family)
        ped = find_pedigree(project, family, project_family)
        rewrite_config_yaml(config_path, project_family=project_family, hpo=hpo, ped=ped)

        # Create samples.tsv / units.tsv
        (family_dir / "samples.tsv").write_text("sample\tBAM\tcase_or_control\n")

        units_tsv = family_dir / "units.tsv"
        units_tsv.write_text("family\tplatform\tsmall_variant_vcf\tpbsv_vcf\tcnv_dir\n")
        with units_tsv.open("a") as out:
            out.write(f"{project_family}\tPACBIO\t{deepvariant}\t{sv}\tcnv/vcfs\n")

        ensure_dir(family_dir / "cnv")
    return deepvariant, sv


def add_sample_inputs(
    *,
    family_dir: Path,
    project: str,
    sequence_id: str,
    project_id_norm: str,
    deepvariant: Path,
    sv: Path,
) -> None:
    # samples.tsv
    bam = FILES_FROM_IRODS / project / f"{sequence_id}.GRCh38.haplotagged.bam"
    if not bam.exists():
        try:
            bam = glob.glob(f"{FILES_FROM_IRODS}/{project}/{sequence_id}*GRCh38.aligned.haplotagged.bam")[0] # older pipeline runs
        except IndexError:
            raise FileNotFoundError(f"No BAM found for {sequence_id}")
    project_sample = project_sample_from_project_id(project_id_norm)
    with (family_dir / "samples.tsv").open("a") as out:
        out.write(f"{project_sample}\t{bam}\n")

    # CNV copy
    cnv_dir = family_dir / "cnv" / "vcfs"
    ensure_dir(cnv_dir)
    cnv_src_exact = FILES_FROM_IRODS / project / f"{sequence_id}.GRCh38.hificnv.vcf.gz"
    if cnv_src_exact.exists():
        copy_with_sidecars(cnv_src_exact, cnv_dir)
    else:
        patt = f"hificnv*{sequence_id}*.vcf.gz"
        matches = sorted((FILES_FROM_IRODS / project).glob(patt))
        if not matches:
            raise FileNotFoundError(f"No CNV VCFs found for {sequence_id} (pattern {patt})")
        for m in matches:
            copy_with_sidecars(m, cnv_dir)

    # Replace sample IDs in VCFs
    # mapping 1: sequence_id -> normalized project_id (deepvariant/sv)
    with tempfile.NamedTemporaryFile("w", delete=False, dir=str(family_dir), prefix="sample_rename_", suffix=".txt") as tf:
        tf.write(f"{sequence_id} {project_id_norm}\n")
        map_seq_to_proj = Path(tf.name)
    try:
        bcftools_reheader_inplace(deepvariant, map_seq_to_proj)
        bcftools_reheader_inplace(sv, map_seq_to_proj)
    finally:
        map_seq_to_proj.unlink()

    # CNV vcfs: map each VCF's existing sample ID -> normalized project_id
    cnv_vcf = Path(glob.glob(f"{cnv_dir}/*{sequence_id}*.vcf.gz")[0])
    cnv_sample = bcftools_query_sample(cnv_vcf)
    with tempfile.NamedTemporaryFile("w", delete=False, dir=str(family_dir), prefix="sample_rename_", suffix=".txt") as tf:
        tf.write(f"{cnv_sample} {project_id_norm}\n")
        map_cnv_to_proj = Path(tf.name)
    try:
        bcftools_reheader_inplace(cnv_vcf, map_cnv_to_proj)
        _run(["tabix", str(cnv_vcf)])
    finally:
        map_cnv_to_proj.unlink()


def validate_pedigrees(analysis_rows: list[AnalysisRow], analyses_path: Path, project: str) -> None:
    # project_family list: third column, cut '.' -f1, uniq
    seen: set[str] = set()
    for r in analysis_rows:
        if r.family_is_header:
            continue
        project_family = r.project_id_raw.split(".", 1)[0]
        if project_family in ("Decoder_ID", "TG_ID"):
            continue
        if not project_family or project_family in seen:
            continue
        seen.add(project_family)

        if ("DSK" in project_family) or ("GYM" in project_family):
            ped_family = project_family.replace("_", "")
        else:
            ped_family = project_family

        ped_matches = sorted((PED_DIR / project).glob(f"{ped_family}*"), key=lambda p: p.stat().st_mtime, reverse=True)
        if not ped_matches:
            raise FileNotFoundError(f"No pedigree found for {project_family} (looked for {ped_family}*)")
        ped = ped_matches[0]

        ped_samples: set[str] = set()
        with ped.open() as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    ped_samples.add(parts[1])

        # Samples in analysis TSV for this family: grep project_family | grep LRWGS | awk '{print $3}' | tr -d '_' | tr '.' '_'
        for rr in analysis_rows:
            if rr.family_is_header:
                continue
            if (project_family not in rr.project_id_raw) or (rr.sample_type != "LRWGS"):
                continue
            sample = _strip_cr(rr.project_id_raw).replace("_", "").replace(".", "_")
            if sample not in ped_samples:
                raise SystemExit(f"Sample {sample} not found in pedigree for {project_family}, exiting")


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="Set up crg2-pacbio analysis directories for PacBio inputs.")
    ap.add_argument("analyses", type=Path, help="Path to sample metadata TSV")
    ap.add_argument("project", help="Project ID, e.g. DECODER")
    ap.add_argument("--creds", default=DEFAULT_CREDS, help="Phenotips credentials CSV (default: PT_credentials.csv)")
    ap.add_argument("--crg2-pacbio", dest="crg2_pacbio", type=Path, default=REPO_ROOT_DEFAULT_CRG2_PACBIO, help="Path to crg2-pacbio repo (default: ~/crg2-pacbio)")
    ap.add_argument("--today", default=None, help="Override date stamp (YYYY-MM-DD). Default: today.")
    args = ap.parse_args(argv)

    if not args.analyses.is_file():
        raise SystemExit(f"Error: sample file does not exist: {args.analyses}")

    _which_or_die("bcftools")
    _which_or_die("tabix")

    today = args.today or _dt.date.today().isoformat()

    analysis_dir = ANALYSES_BASE / args.project
    ensure_dir(analysis_dir)

    rows = parse_analysis_tsv(args.analyses)

    # Download HPO + pedigrees from Phenotips (mirrors bash)
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
        if r.family_is_header:
            continue

        family = r.family
        sequence_id = _strip_cr(r.sequence_id)

        try:
            project_id_norm = normalize_project_id(family, r.project_id_raw)
        except ValueError as e:
            raise SystemExit(str(e)) from e

        project_family = project_family_from_project_id(project_id_norm)
        family_dir = analysis_dir / project_family / f"PacBio_{today}"

        deepvariant, sv = setup_family_once(
            family=family,
            sequence_id=sequence_id,
            project=args.project,
            project_id_norm=project_id_norm,
            project_family=project_family,
            family_dir=family_dir,
            crg2_pacbio=args.crg2_pacbio,
            today=today,
        )

        add_sample_inputs(
            family_dir=family_dir,
            project=args.project,
            sequence_id=sequence_id,
            project_id_norm=project_id_norm,
            deepvariant=deepvariant,
            sv=sv,
        )

    validate_pedigrees(rows, args.analyses, args.project)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

