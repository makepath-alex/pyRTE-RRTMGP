#!/usr/bin/env python3
import json
import numpy as np
import pyrte_rrtmgp.pyrte_rrtmgp as py
import os
import copy as cp
import zipfile

def load_data(path):

    # Read dimensions
    with open(os.path.join(path, "inputs-dim.bin"), "rb") as f:
        dims = np.fromfile(f, dtype=np.int32, count=9)
    ncol, nlay, nbnd, ngpt, nflav, neta, npres, ntemp, nPlanckTemp = dims

    print("Dimensions:", ncol, nlay, nbnd, ngpt, nflav, neta, npres, ntemp, nPlanckTemp)

    # Read integer data
    data = {
        "ncol": ncol,
        "nlay": nlay,
        "nbnd": nbnd,
        "ngpt": ngpt,
        "nflav": nflav,
        "neta": neta,
        "npres": npres,
        "ntemp": ntemp,
        "nPlanckTemp": nPlanckTemp,
        "tlay": np.fromfile(os.path.join(path, "inputs-tlay.bin"), dtype=np.float64, count=ncol * nlay),#.reshape((ncol, nlay), order='C'),
        "tlev": np.fromfile(os.path.join(path, "inputs-tlev.bin"), dtype=np.float64, count=ncol * (nlay + 1)),#.reshape((ncol, nlay + 1), order='C'),
        "tsfc": np.fromfile(os.path.join(path, "inputs-tsfc.bin"), dtype=np.float64, count=ncol),
        "sfc_lay": np.fromfile(os.path.join(path, "inputs-sfc_lay.bin"), dtype=np.int32, count=1)[0],
        "fmajor": np.fromfile(os.path.join(path, "inputs-fmajor.bin"), dtype=np.float64, count=2 * 2 * 2 * ncol * nlay * nflav),#.reshape((2, 2, 2, ncol, nlay, nflav), order='C'),
        "jeta": np.fromfile(os.path.join(path, "inputs-jeta.bin"), dtype=np.int32, count=2 * ncol * nlay * nflav),#.reshape((2, ncol, nlay, nflav), order='C'),
        "tropo": np.fromfile(os.path.join(path, "inputs-tropo.bin"), dtype=np.bool_, count=ncol * nlay),#.reshape((ncol, nlay), order='C'),
        "jtemp": np.fromfile(os.path.join(path, "inputs-jtemp.bin"), dtype=np.int32, count=ncol * nlay),#.reshape((ncol, nlay), order='C'),
        "jpress": np.fromfile(os.path.join(path, "inputs-jpress.bin"), dtype=np.int32, count=ncol * nlay),#.reshape((ncol, nlay), order='C'),
        "gpoint_bands": np.fromfile(os.path.join(path, "inputs-gpoint_bands.bin"), dtype=np.int32, count=ngpt),
        "band_lims_gpt": np.fromfile(os.path.join(path, "inputs-band_lims_gpt.bin"), dtype=np.int32, count=2 * nbnd),#.reshape((2, nbnd), order='C'),
        "pfracin": np.fromfile(os.path.join(path, "inputs-pfracin.bin"), dtype=np.float64, count=ntemp * neta * (npres + 1) * ngpt),#.reshape((ntemp, neta, npres + 1, ngpt), order='C'),
        "temp_ref_min": np.fromfile(os.path.join(path, "inputs-temp_ref_min.bin"), dtype=np.float64, count=1)[0],
        "totplnk_delta": np.fromfile(os.path.join(path, "inputs-totplnk_delta.bin"), dtype=np.float64, count=1)[0],
        "totplnk": np.fromfile(os.path.join(path, "inputs-totplnk.bin"), dtype=np.float64, count=nPlanckTemp * nbnd),#.reshape((nPlanckTemp, nbnd), order='C'),
        "gpoint_flavor": np.fromfile(os.path.join(path, "inputs-gpoint_flavor.bin"), dtype=np.int32, count=2 * ngpt),#.reshape((2, ngpt), order='C'),
        "sfc_src": np.fromfile(os.path.join(path, "inputs-sfc_src.bin"), dtype=np.float64, count=ncol * ngpt),#.reshape((ncol, ngpt), order='C'),
        "lay_src": np.fromfile(os.path.join(path, "inputs-lay_src.bin"), dtype=np.float64, count=ncol * nlay * ngpt),#.reshape((ncol, nlay, ngpt), order='C'),
        "lev_src": np.fromfile(os.path.join(path, "inputs-lev_src.bin"), dtype=np.float64, count=ncol * (nlay + 1) * ngpt),#.reshape((ncol, nlay + 1, ngpt), order='C'),
        "sfc_source_Jac": np.fromfile(os.path.join(path, "inputs-sfc_source_Jac.bin"), dtype=np.float64, count=ncol * ngpt),#.reshape((ncol, ngpt), order='C')
    }

    return data

def test_rrtmgp_compute_Planck_source():
    """
    Test the Fortran function `rrtmgp_compute_Planck_source` by comparing its outputs
    against a reference.
    """
    output_keys = [
        "sfc_src",
        "lay_src",
        "lev_src",
        "sfc_source_Jac"
    ]

    path = "/data/vkm/code/makepath/rte-rrtmgp/build/examples/rfmip-clear-sky/"
    pathRef = "/data/vkm/code/makepath/rte-rrtmgp/tmp"

    # Load input data
    input_data  = load_data(path)
    output_data = load_data(pathRef)

    # Prepare arguments for Fortran function
    args = list(input_data.values())

    # Call the Fortran function
    py.rrtmgp_compute_Planck_source(*args)

    # Validate outputs
    for key in output_keys:
        computed  = input_data[key]
        reference = output_data[key]

        threshold = 0.0#7.e-4

        diff  = np.abs(reference - computed)
        mask  = diff > threshold
        count = np.sum(mask)

        if count > 0:
            indices = np.argwhere(mask)
            for idx in indices:
                ref_val = reference[tuple(idx)]
                computed_val = computed[tuple(idx)]
                diff_val = diff[tuple(idx)]
                # print(f"#{idx}\t{ref_val}\t| {computed_val}\t>> {diff_val}")

        print(f"{key} shape: {computed.shape}\tDiff count: {count}\tMax diff: {np.max(diff)}")

        # assert np.allclose(ref, computed, 0.0), f"Mismatch in {key}: {np.max(np.abs(ref - computed))} difference"

test_rrtmgp_compute_Planck_source()
