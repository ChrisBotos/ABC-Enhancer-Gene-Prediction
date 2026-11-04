import os
import subprocess
import tempfile
import unittest

import pandas as pd

CONFIG_FILE = "tests/config/generic_config.yml"
DHS_BAM = "example_chr/chr22/ENCFF860XAE.chr22.sorted.se.bam"
H3K27AC_BAM = "example_chr/chr22/ENCFF790GFL.chr22.sorted.se.bam"


def dry_run_commands(biosample_row: dict, outdir: str) -> str:
    """
    Writes a one-row biosamples table and returns the shell commands of a dry run.
    """
    biosamples_file = os.path.join(outdir, "biosamples.tsv")
    row = {
        "biosample": "K562_chr22",
        "DHS": DHS_BAM,
        "ATAC": "",
        "H3K27ac": H3K27AC_BAM,
        "default_accessibility_feature": "DHS",
        "HiC_file": "",
        "HiC_type": "",
        "HiC_resolution": "",
        "alt_TSS": "",
        "alt_genes": "",
        **biosample_row,
    }
    pd.DataFrame([row]).to_csv(biosamples_file, sep="\t", index=False)
    cmd = [
        "snakemake",
        "-n",
        "-p",
        "-F",
        "--configfile",
        CONFIG_FILE,
        "--config",
        f"biosamplesTable={biosamples_file}",
        f"results_dir={outdir}/results/",
    ]
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def macs2_command(commands: str) -> str:
    start = commands.index("macs2 callpeak")
    return commands[start : commands.index("\n\n", start)]


class TestElementSource(unittest.TestCase):
    def test_default_calls_peaks_on_accessibility(self) -> None:
        with tempfile.TemporaryDirectory() as outdir:
            commands = dry_run_commands({}, outdir)
        macs2 = macs2_command(commands)
        self.assertIn(f"-t {DHS_BAM}", macs2)
        self.assertIn('CONTROL=""', commands)
        self.assertIn(f"--accessibility {DHS_BAM}", commands)

    def test_element_source_replaces_accessibility_for_candidates(self) -> None:
        row = {"element_source": H3K27AC_BAM, "element_source_control": DHS_BAM}
        with tempfile.TemporaryDirectory() as outdir:
            commands = dry_run_commands(row, outdir)
        macs2 = macs2_command(commands)
        self.assertIn(f"-t {H3K27AC_BAM}", macs2)
        self.assertIn(f'CONTROL="-c {DHS_BAM}"', commands)
        self.assertIn(f"--accessibility {H3K27AC_BAM}", commands)
        # Neighborhoods still count activity on the DHS column.
        self.assertIn(f"--DHS {DHS_BAM}", commands)


if __name__ == "__main__":
    unittest.main()
