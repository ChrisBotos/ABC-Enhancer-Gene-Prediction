import os
import sys
import tempfile
import unittest

import pandas as pd

SCRIPTS_DIR = os.path.abspath("workflow/scripts")
sys.path.insert(0, SCRIPTS_DIR)
from neighborhoods import annotate_genes_with_features, load_genes

EXPECTED_DIR = "tests/expected_output/generic/K562_chr22"
GENES = os.path.join(EXPECTED_DIR, "processed_genes_file.bed")
CHROM_SIZES = "reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv"
UBIQUITOUS_GENES = "reference/UbiquitouslyExpressedGenes.txt"
DHS_BAM = "example_chr/chr22/ENCFF860XAE.chr22.sorted.se.bam"
H3K27AC_BAM = "example_chr/chr22/ENCFF790GFL.chr22.sorted.se.bam"


def gene_list(outdir: str, features: dict, expression_features: dict) -> pd.DataFrame:
    chrom_sizes_bed = os.path.join(outdir, "chrom_sizes.bed")
    sizes = pd.read_csv(CHROM_SIZES, sep="\t", header=None, names=["chr", "end"])
    sizes.assign(start=0)[["chr", "start", "end"]].to_csv(
        chrom_sizes_bed, sep="\t", header=False, index=False
    )
    genes, _ = load_genes(
        file=GENES,
        ue_file=UBIQUITOUS_GENES,
        chrom_sizes=CHROM_SIZES,
        outdir=outdir,
        expression_table_list="",
        gene_id_names="symbol",
        primary_id="symbol",
        cellType=None,
        class_gene_file=None,
    )
    annotate_genes_with_features(
        genes=genes,
        genome_sizes=CHROM_SIZES,
        genome_sizes_bed=chrom_sizes_bed,
        chrom_sizes_map=sizes.set_index("chr")["end"].to_dict(),
        features=features,
        outdir=outdir,
        default_accessibility_feature="DHS",
        expression_features=expression_features,
    )
    return pd.read_csv(os.path.join(outdir, "GeneList.txt"), sep="\t")


class TestPromoterActivityQuantile(unittest.TestCase):
    def setUp(self) -> None:
        self.reference = {"DHS": [DHS_BAM], "H3K27ac": [H3K27AC_BAM]}

    def test_default_matches_expected_gene_list(self) -> None:
        with tempfile.TemporaryDirectory() as outdir:
            test = gene_list(outdir, self.reference, {})
        expected = pd.read_csv(os.path.join(EXPECTED_DIR, "Neighborhoods/GeneList.txt"), sep="\t")
        pd.testing.assert_series_equal(
            test["PromoterActivityQuantile"], expected["PromoterActivityQuantile"]
        )

    def test_reference_expression_tracks_leave_gene_list_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as outdir:
            default = gene_list(outdir, self.reference, {})
        with tempfile.TemporaryDirectory() as outdir:
            filled = gene_list(outdir, self.reference, self.reference)
        pd.testing.assert_frame_equal(filled, default)

    def test_expression_tracks_hold_the_call_when_a_slot_is_replaced(self) -> None:
        replaced = {"DHS": [H3K27AC_BAM], "H3K27ac": [DHS_BAM]}
        with tempfile.TemporaryDirectory() as outdir:
            default = gene_list(outdir, self.reference, {})
        with tempfile.TemporaryDirectory() as outdir:
            swapped = gene_list(outdir, replaced, {})
        with tempfile.TemporaryDirectory() as outdir:
            held = gene_list(outdir, replaced, self.reference)
        pd.testing.assert_series_equal(
            held["PromoterActivityQuantile"], default["PromoterActivityQuantile"]
        )
        self.assertFalse(
            swapped["H3K27ac.RPM.TSS1Kb"].equals(default["H3K27ac.RPM.TSS1Kb"])
        )
        pd.testing.assert_series_equal(
            held["H3K27ac.RPM.TSS1Kb"], swapped["H3K27ac.RPM.TSS1Kb"]
        )


if __name__ == "__main__":
    unittest.main()
