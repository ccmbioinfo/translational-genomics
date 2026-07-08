
# Archiving data to iRODS

`01_irsync_to_iRODS.sh`

This script archives directories on Isilon/PowerScale to iRODS using `irsync` with the `-K` flag to perform checksums. Likely overkill, but it performs `irsync` twice to verify that the files have been archived appropriately. Replace `your.email@email.ca` in the SBATCH directive with your email to receive an email upon completion or failure of the job.

Usage: ` 01_irsync_to_iRODS.sh <local_path> <archival_path>`

`02_final_check.sh`

This script compares the numbers of files in the directory on Isilon/PowerScale to the corresponding directory on iRODS. The number of files in the local directory, multiplied by two, should equal the number of files in the iRODS directory (each file in iRODS is mirrored). If the check is `OK`, you can proceed with removing the files locally that you have just archived. If the file outputs a `WARNING`, try rsyncing your files again.

Usage: `02_final_check.sh <local_path> <archival_path>`

Example output:

| iRODS_count_with_mirror |	Isilon_count_x2 | check |
| ----------------------- | --------------- | ------|
| 754 | 754 | OK | 