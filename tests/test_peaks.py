import os
import sys
import tempfile
import unittest

import pandas as pd

SCRIPTS_DIR = os.path.abspath("workflow/scripts")
sys.path.insert(0, SCRIPTS_DIR)
from peaks import get_peak_enrichment, make_candidate_regions_from_summits

EXPECTED_PEAKS_DIR = "tests/expected_output/generic/K562_chr22/Peaks"
NARROW_PEAK = os.path.join(EXPECTED_PEAKS_DIR, "macs2_peaks.narrowPeak.sorted")
CHROM_SIZES = "reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv"
BLOCKLIST = "reference/hg38/GRCh38_unified_blacklist.bed"
INCLUDELIST = "example_chr/chr22/RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.hg38.TSS500bp.bed"
DHS_BAM = "example_chr/chr22/ENCFF860XAE.chr22.sorted.se.bam"


def make_candidates(outdir: str, rank_by: str, n_enhancers: int) -> str:
    chrom_sizes_bed = os.path.join(outdir, "chrom_sizes.bed")
    sizes = pd.read_csv(CHROM_SIZES, sep="\t", header=None, names=["chr", "end"])
    sizes.assign(start=0)[["chr", "start", "end"]].to_csv(
        chrom_sizes_bed, sep="\t", header=False, index=False
    )
    make_candidate_regions_from_summits(
        macs_peaks=NARROW_PEAK,
        accessibility_files=[DHS_BAM],
        genome_sizes=CHROM_SIZES,
        genome_sizes_bed=chrom_sizes_bed,
        regions_includelist=INCLUDELIST,
        regions_blocklist=BLOCKLIST,
        n_enhancers=n_enhancers,
        peak_extend=250,
        outdir=outdir,
        rank_by=rank_by,
    )
    return os.path.join(outdir, os.path.basename(NARROW_PEAK) + ".candidateRegions.bed")


class TestPeaks(unittest.TestCase):
    def test_rank_by_reads_matches_expected_candidates(self) -> None:
        with tempfile.TemporaryDirectory() as outdir:
            with open(make_candidates(outdir, "reads", 150000)) as f:
                test_contents = f.read()
        with open(NARROW_PEAK + ".candidateRegions.bed") as f:
            expected_contents = f.read()
        self.assertEqual(test_contents, expected_contents)

    def test_enrichment_scores_are_narrowpeak_fold_enrichment(self) -> None:
        with tempfile.TemporaryDirectory() as outdir:
            scores = pd.read_csv(
                get_peak_enrichment(NARROW_PEAK, outdir), sep="\t", header=None
            )
        peaks = pd.read_csv(NARROW_PEAK, sep="\t", header=None)
        pd.testing.assert_frame_equal(
            scores, peaks[[0, 1, 2, 6]].set_axis([0, 1, 2, 3], axis=1)
        )

    def test_rank_by_enrichment_keeps_the_most_enriched_peak(self) -> None:
        peaks = pd.read_csv(NARROW_PEAK, sep="\t", header=None)
        top = peaks.loc[peaks[6].idxmax()]
        summit = top[1] + top[9]
        with tempfile.TemporaryDirectory() as outdir:
            candidates = pd.read_csv(
                make_candidates(outdir, "enrichment", 1), sep="\t", header=None
            )
        self.assertTrue(
            ((candidates[1] <= summit) & (candidates[2] >= summit)).any()
        )


if __name__ == "__main__":
    unittest.main()
