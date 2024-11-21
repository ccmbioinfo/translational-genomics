from datetime import date
import sys

# parse through iRods archive files/directories (i.e. 'ils -r /resarchivezone/tcag/translational_genomics/pacbio/decoder/ > iRods_<date>.txt')
# get PacBio deepvariant VCF, pbsv VCF, BAMs, and CNV reports
# usage: python3 parse_irods_pacbio.py iRods_<date>.txt family1,family2,family3

irods = sys.argv[1]
families = sys.argv[2]
families = list(families.split(","))
deepvariant_dict = {}
pbsv_dict = {}
bam_dict = {}
trgt_dict = {}
cnv_dict = {}
with open(irods, 'r') as file:
    for line in file: 
        line = line.strip()
        for family in families:
            family_DNA = family.split("_")[0] + "_DNA_" + family.split("_")[1]
            if family in line or family_DNA in line:
                if 'fam' in line and 'humanwgs' in line and ':' in line:
                    family_dir = line.split('/')[-2]
                    family_prefix = family_dir.replace(".joint", "")
                    family = line.split('/')[-2].split('-')[0]
                    irods_dir = line.replace(":", "")
                    pbsv = f"{irods_dir}/{family_prefix}.joint.GRCh38.pbsv.phased.vcf.gz"
                    pbsv_dict[family] = pbsv
                    deepvariant =  f"{irods_dir}/{family_prefix}.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz"
                    deepvariant_dict[family] = deepvariant
                    print(family, pbsv, deepvariant)
                elif 'GRCh38.aligned.haplotagged.bam' in line and 'bai' not in line:
                    sample = line.split('.')[0]
                    bam = f"{irods_dir}/{line}"
                    print(bam)
                    bam_dict[sample] = bam
                elif '.trgt.sorted.vcf.gz' in line:
                    sample = line.split('.')[0]
                    trgt = f"{irods_dir}/{line}"
                    print(trgt)
                    trgt_dict[sample] = trgt
                elif 'annotations:' in line and 'CPG' not in line and 'fam' not in line:
                    sample = line.split('/')[-2]
                    cnv_irods_dir = line.replace(":", "")
                    cnv = f"{cnv_irods_dir}/{sample}.hificnv_annPipelineRev1.6.0_20230818_hg38_PACBIO_cnv.tagged.tsv"
                    print(cnv)
                    cnv_dict[sample] = cnv

today = date.today().strftime("%Y-%m-%d")
with open(f'../sample_sheets/deepvariant_{today}.txt', 'w') as file:
    for family, vcf in deepvariant_dict.items():
        file.write(f"{vcf}\n")

with open(f'../sample_sheets/pbsv_{today}.txt', 'w') as file:
    for family, vcf in pbsv_dict.items():
        file.write(f"{vcf}\n")

with open(f'../sample_sheets/BAMs_{today}.txt', 'w') as file:
    for sample, bam in bam_dict.items():
        file.write(f"{bam}\n")
        file.write(f"{bam}.bai\n")
        trgt_bam = bam.replace(".bam", "") + ".trgt.spanning.sorted.bam"
        file.write(f"{trgt_bam}\n")
        file.write(f"{trgt_bam}.bai\n")

with open(f'../sample_sheets/trgt_{today}.txt', 'w') as file:
    for sample, trgt in trgt_dict.items():
        file.write(f"{trgt}\n")

with open(f'../sample_sheets/CNV_{today}.txt', 'w') as file:
    for sample, cnv in cnv_dict.items():
        file.write(f"{cnv}\n")
                


