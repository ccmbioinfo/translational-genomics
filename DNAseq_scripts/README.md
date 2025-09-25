## DNAseq Scripts Documentation

## The scripts in this folder perform a variety of tasks related to processing exomes, short-read genomes, and long-read genomes for the SickKids Translational Genomics Node. 

### copy_files_from_irods.sh
A bash script that copies files from the iRods archive to Isilon storage. Takes two arguments:
- `irods_files`: A file containing list of files to copy from iRods
- `project`: Project name/directory to copy files to (e.g. DECODER)

Usage:
`sh copy_files_from_irods.sh <irods_files> <project>`

### get_HPO_pedigree_genome_clinic.py 
A Python script that retrieves pedigree information and HPO terms from Genome Clinic Phenotips. Contains functions to:
- Get pedigree information including family relationships
- Map pedigree nodes to C4R IDs
- Retrieve affected status and sex information
- Get parental relationships
- Takes two arguments: 
    - `-sample_sheet`: Tab-separated sample sheet with at minimum Family_ID, Sample_ID, and Decoder_ID columns
    - `-credentials`: CSV file containing Phenotips username and password

Usage:
`python3 get_HPO_pedigree_genome_clinic.py  -sample_sheet <sample sheet TSV> -credentials <credentials CSV>`

### HPO_excel_to_text.py
Converts HPO terms stored in Excel format to text files. Specifically:
- Reads HPO terms from Excel sheets for clinical samples
- Cleans data by removing header rows and empty entries
- Outputs tab-separated text files with gene symbols, IDs, and HPO terms

Usage:
`python3 HPO_excel_to_text.py`

### organize_crg2-pacbio_reports.sh
Organizes analysis reports from the crg2-pacbio pipeline:
- Copies and organizes various report types (coding, panel, CNV, SV, etc.)
- Adds HPO terms to small variant reports
- Consolidates reports into a single directory
- One argument: `analysis path`, e.g. DSK002/PacBio

Usage:
`sh organize_crg2-pacbio_reports.sh <analysis_path>`

### query_files_from_irods.sh
A bash script that queries the iRods paths for crg2-pacbio input files, such as the small variant VCFs, BAMs, etc. 
- Takes two arguments:
    - `LIMSID`: sequencing batch ID, e.g. COS30304
    - `project`: i.e. DECODER or genesteps
- Outputs filepaths to /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${PROJECT}/iRods_paths/${LIMSID}_PacBio_files.txt

Usage: 
`sh query_files_from_irods.sh <LIMSID> <project>`

### PacBio_setup.sh
A bash script that sets up analysis directories and configuration for PacBio data processing:
- Takes three arguments:
    - `analyses`: Tab-separated sample sheet with Family_ID, Sample_ID, and Decoder_ID columns
    - `iRods_query_date`: Date when iRods files were queried (YYYY-MM-DD format)
    - `project`: Project name (e.g., genesteps or DECODER)
- Sets up analysis directories with required pipeline files
- Configures HPO terms and pedigree information
- Handles sample renaming and file organization

Usage:
`sh PacBio_setup.sh <analyses TSV> <iRods_date> <project>`

### parse_irods_pacbio.py
Parses iRods archive files for PacBio DECODER project data:
- Takes two arguments:
    - iRods file listing (from `ils -r` command)
    - Comma-separated list of family IDs
- Extracts paths for:
    - DeepVariant VCFs
    - PBSV VCFs
    - BAM files
    - TRGT files
    - CNV reports
- Outputs dated file lists for each data type

Usage:
`python3 parse_irods_pacbio.py <iRods_file> <family1,family2,family3>`

### parse_irods_pacbio_genesteps.py
Similar to parse_irods_pacbio.py but specifically for Gene-STEPS project data:
- Takes one argument: iRods file listing
- Automatically processes all families in the Gene-STEPS iRods directory
- Creates dated file lists for all data types

Usage:
`python3 parse_irods_pacbio_genesteps.py <iRods_file>`

### run_TRGT_repeat_outliers_and_denovo.sh
Sets up and configures TRGT outlier and de novo repeat analyses for a family:
- Takes two arguments:
    - `family`: Family ID
    - `project`: Project name
- Configures HPO terms and pedigree information
- Sets up TRGT-denovo analysis

Usage:
`sh run_TRGT_repeat_outliers_and_denovo.sh <family> <project>`

### visualize_repeat.sh
Creates visualization plots for repeat regions:
- Takes three arguments:
    - `VCF`: Input VCF file
    - `BAM`: BAM file with spanning reads
    - `REPEAT_ID`: ID of repeat to visualize
- Generates SVG plots showing repeat structure and methylation
- Uses TRGT plotting functionality

Usage:
`sh visualize_repeat.sh <VCF> <BAM> <REPEAT_ID>`

### WGS_reanalysis_DPLM_input.sh
Sets up CRG2 pipeline reanalysis for sequences (WES or WGS) from SickKids DPLM clinical genetics lab:
- Takes sample TSV file and project name (e.g. DECODER) as input
- Validates input files and directories
- Processes sample IDs and sets up analysis directories
- Handles DECODER ID formatting for pipeline compatibility

Usage:
`sh WGS_reanalysis_DPLM_input.sh <sample TSV> <project name>`