from typing import TYPE_CHECKING

from STMiner.IO.read_bmk import read_bmk
from STMiner.IO.read_h5ad import read_h5ad
from STMiner.IO.read_stereo import read_gem_file


if TYPE_CHECKING:
    from STMiner.SPFinder import SPFinder


def read_h5ad_for_spfinder(
    spfinder: "SPFinder",
    file,
    amplification=1,
    bin_size=1,
    merge_bin=False,
):
    spfinder.set_adata(
        read_h5ad(
            file,
            amplification=amplification,
            bin_size=bin_size,
            merge_bin=merge_bin,
        )
    )
    spfinder._scope = (0, max(spfinder.adata.obs["y"].max(), spfinder.adata.obs["x"].max()))


def read_bmk_dir_for_spfinder(spfinder: "SPFinder", file, bin_size=1):
    spfinder.set_adata(read_bmk(file, bin_size=bin_size))


def read_gem_for_spfinder(spfinder: "SPFinder", file, bin_size=40):
    spfinder.set_adata(read_gem_file(file, bin_size=bin_size))
