import glob
import os
import subprocess

# Find all CSV files
for csv in glob.glob("../Data/*.csv"):
    basename = os.path.basename(csv)  
    basename_no_prefix = basename.replace("PCs_", "")
    arg1 = basename_no_prefix.split("_")[0]
    arg2 = os.path.splitext(basename_no_prefix)[0]

    # Count columns in header
    with open(csv) as f:
        header = f.readline().strip()
    n_col = len(header.split(","))

    # Loop over columns
    for col in range(3, n_col + 1):
        pc = col - 2
        cmd = [
            "sbatch",
            "./master_script.sh",
            arg1,
            f"{arg2}_PC{pc}",
            csv,
            str(col),
        ]
        subprocess.run(cmd)

