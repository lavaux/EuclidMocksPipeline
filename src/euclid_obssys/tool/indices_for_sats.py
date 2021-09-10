################################################################################
### Author: Tiago Castro, Pierluigi Monaco, Guilhem Lavaux                   ###
###                                                                          ###
################################################################################
from . import register_tool
import sys
from ..config import readConfig


@register_tool
def createIndicesForSats(config: str) -> None:
    from astropy.io import fits
    from astropy.table import Table
    import numpy as np
    from euclid_obssys.config import readConfig
    import sys
    import zarr

    input = readConfig(config)

    print("# Running createIndicesForSats.py")

    source = zarr.open_group(input.master_fname(), mode="r")
    cat = source["catalog"]

    # cat=fits.getdata('../Products/RawCatalogs/8614_100sqdeg.fits')
    # cat=fits.getdata('../Products/catalog.fits')

    print(f"read {input.master_fname()}")

    halo_id = cat["halo_id"]
    kind = cat["kind"]
    index = np.argsort(halo_id)

    print("sorting done")

    Ngal = len(cat)
    halo_index = np.zeros_like(halo_id)
    i = 0
    wmh = 0
    step = 1e6
    while i < Ngal:

        if i > step:
            print("progress: %d out of %d" % (i, Ngal))
            step += 1e6

        this_halo = halo_id[index[i]]
        first = i
        last = i + 1

        if last == Ngal:
            break

        while halo_id[index[last]] == this_halo:
            last += 1
            if last == Ngal:
                break

        # print("i=%d, this_halo=%d, first=%d, last=%d"%(i,this_halo,first,last))

        if first > last + 1:
            print("gal %d, section: %d-%d" % (i, first, last))

        if last == i + 1:
            halo_index[index[i]] = index[i]
            # print("only one")

        else:
            section = kind[index[first:last]]
            main = np.where(section == 0)[0]
            if i < 1000:
                print(
                    "main: ",
                    main,
                    first,
                    last,
                    kind[index[first]],
                    kind[index[first + main[0]]],
                )

            if len(main) > 0:
                halo_index[index[first:last]] = index[first + main[0]]
            else:
                halo_index[index[first:last]] = -1
                wmh += 1
                # print("section without main halo")

        i = last

    print("loop done")

    # t=Table([halo_index], names=('indices',))
    # fname='../Products/RawCatalogs/8614_100sqdeg_indices.fits'
    # t.write(fname,format='fits',overwrite=True)

    out_index = zarr.open_group(input.indices_fname())
    out_index["indices"] = halo_index

    print("written file " + input.indices_fname())
