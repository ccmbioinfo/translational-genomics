import datetime
import glob
import logging
import pandas as pd
import os
import pathlib

data_dir = "/hpf/largeprojects/tgnode/data"
trash_dir = "/hpf/largeprojects/tgnode/trash"
sample_sheet = "/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/sample_sheets/DECODER/DECODER_analyses_all.tsv"
log_file = f"/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/logs/cleanup_{datetime.datetime.now().strftime('%Y-%m-%d')}.log"

logging.basicConfig(filename=log_file, level=logging.INFO)

# Read sample sheet
sample_sheet = pd.read_csv(sample_sheet, sep="\t")
sample_sheet_to_delete = sample_sheet[(sample_sheet["Status"] == "Done") & (sample_sheet["Lab"].isin(["DPLM", "GeneDx", "Prevention_Genetics"]))]

# Get list of sequence IDs
for _, row in sample_sheet_to_delete.iterrows():
    sequence_ID = row["Sample_ID"]
    DECODER_ID = row["Decoder_ID"]
    files = glob.glob(f"{data_dir}/**/{sequence_ID}*", recursive=True)
    if len(files) == 0:
        logging.info(f"No files found for {DECODER_ID} {sequence_ID}")
        continue
    for file in files:
        if os.path.exists(file):
            logging.info(f"{DECODER_ID}: moving {file} to {trash_dir}")
            os.rename(file, f"{trash_dir}/{os.path.basename(file)}")
        else:
            logging.info(f"File {file} does not exist")