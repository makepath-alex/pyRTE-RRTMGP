# #!/usr/bin/env python3
# import json
# import numpy as np
# import pyrte_rrtmgp.pyrte_rrtmgp as py
# import os
# import copy as cp

# def load_data(path):

#     # Read dimensions
#     with open(os.path.join(path, "inputs-dims.bin"), "rb") as f:
#         dims = np.fromfile(f, dtype=np.int32, count=9)
#     ncol, nlay, nbnd, ngpt, ngas, nflav, neta, npres, ntemp = dims

#     with open(os.path.join(path, "inputs-nminors.bin"), "rb") as f:
#         nminors = np.fromfile(f, dtype=np.int32, count=4)
#     nminorlower, nminorklower, nminorupper, nminorkupper = nminors

#     print("Dimensions:", ncol, nlay, nbnd, ngpt, ngas, nflav, neta, npres, ntemp)
#     print("NMinors:", nminorlower, nminorklower, nminorupper, nminorkupper)

#     # Read integer data
#     data = {
#         "ncol": ncol,
#         "nlay": nlay,
#         "nbnd": nbnd,
#         "ngpt": ngpt,
#         "ngas": ngas,
#         "nflav": nflav,
#         "neta": neta,
#         "npres": npres,
#         "ntemp": ntemp,
#         "nminorlower": nminorlower,
#         "nminorklower": nminorklower,
#         "nminorupper": nminorupper,
#         "nminorkupper": nminorkupper,
#         "idx_h2o": np.fromfile(os.path.join(path, "inputs-idx_h2o.bin"), dtype=np.int32, count=1)[0],
#         "gpoint_flavor": np.fromfile(os.path.join(path, "inputs-gpoint_flavor.bin"), dtype=np.int32, count=2 * ngpt),
#         "band_lims_gpt": np.fromfile(os.path.join(path, "inputs-band_lims_gpt.bin"), dtype=np.int32, count=2 * nbnd),
#         "kmajor": np.fromfile(os.path.join(path, "inputs-kmajor.bin"), dtype=np.float64, count=ntemp * neta * (npres+1) * ngpt),
#         "kminor_lower": np.fromfile(os.path.join(path, "inputs-kminor_lower.bin"), dtype=np.float64, count=ntemp * neta * nminorklower),
#         "kminor_upper": np.fromfile(os.path.join(path, "inputs-kminor_upper.bin"), dtype=np.float64, count=ntemp * neta * nminorkupper),
#         "minor_limits_gpt_lower": np.fromfile(os.path.join(path, "inputs-minor_limits_gpt_lower.bin"), dtype=np.int32, count=2 * nminorlower),
#         "minor_limits_gpt_upper": np.fromfile(os.path.join(path, "inputs-minor_limits_gpt_upper.bin"), dtype=np.int32, count=2 * nminorupper),
#         "minor_scales_with_density_lower": np.fromfile(os.path.join(path, "inputs-minor_scales_with_density_lower.bin"), dtype=np.bool_, count=nminorlower),
#         "minor_scales_with_density_upper": np.fromfile(os.path.join(path, "inputs-minor_scales_with_density_upper.bin"), dtype=np.bool_, count=nminorupper),
#         "scale_by_complement_lower": np.fromfile(os.path.join(path, "inputs-scale_by_complement_lower.bin"), dtype=np.bool_, count=nminorlower),
#         "scale_by_complement_upper": np.fromfile(os.path.join(path, "inputs-scale_by_complement_upper.bin"), dtype=np.bool_, count=nminorupper),
#         "idx_minor_lower": np.fromfile(os.path.join(path, "inputs-idx_minor_lower.bin"), dtype=np.int32, count=nminorlower),
#         "idx_minor_upper": np.fromfile(os.path.join(path, "inputs-idx_minor_upper.bin"), dtype=np.int32, count=nminorupper),
#         "idx_minor_scaling_lower": np.fromfile(os.path.join(path, "inputs-idx_minor_scaling_lower.bin"), dtype=np.int32, count=nminorlower),
#         "idx_minor_scaling_upper": np.fromfile(os.path.join(path, "inputs-idx_minor_scaling_upper.bin"), dtype=np.int32, count=nminorupper),
#         "kminor_start_lower": np.fromfile(os.path.join(path, "inputs-kminor_start_lower.bin"), dtype=np.int32, count=nminorlower),
#         "kminor_start_upper": np.fromfile(os.path.join(path, "inputs-kminor_start_upper.bin"), dtype=np.int32, count=nminorupper),
#         "tropo": np.fromfile(os.path.join(path, "inputs-tropo.bin"), dtype=np.bool_, count=ncol * nlay),
#         "col_mix": np.fromfile(os.path.join(path, "inputs-col_mix.bin"), dtype=np.float64, count=2 * ncol * nlay * nflav),
#         "fmajor": np.fromfile(os.path.join(path, "inputs-fmajor.bin"), dtype=np.float64, count=2 * 2 * 2 * ncol * nlay * nflav),
#         "fminor": np.fromfile(os.path.join(path, "inputs-fminor.bin"), dtype=np.float64, count=2 * 2 * ncol * nlay *nflav),
#         "play": np.fromfile(os.path.join(path, "inputs-play.bin"), dtype=np.float64, count=ncol * nlay),
#         "tlay": np.fromfile(os.path.join(path, "inputs-tlay.bin"), dtype=np.float64, count=ncol * nlay),
#         "col_gas": np.fromfile(os.path.join(path, "inputs-col_gas.bin"), dtype=np.float64, count=ncol * nlay * (ngas + 1)),
#         "jeta": np.fromfile(os.path.join(path, "inputs-jeta.bin"), dtype=np.int32, count=2 * ncol * nlay * nflav),
#         "jtemp": np.fromfile(os.path.join(path, "inputs-jtemp.bin"), dtype=np.int32, count=ncol * nlay),
#         "jpress": np.fromfile(os.path.join(path, "inputs-jpress.bin"), dtype=np.int32, count=ncol * nlay),
#         "tau": np.fromfile(os.path.join(path, "inputs-tau.bin"), dtype=np.float64, count=ncol * nlay *ngpt)
#     }

#     return data

# def load_data1(path):
#     # Read dimensions
#     with open(os.path.join(path, "inputs-dims.bin"), "rb") as f:
#         dims = np.fromfile(f, dtype=np.int32, count=9)
#     ncol, nlay, nbnd, ngpt, ngas, nflav, neta, npres, ntemp = dims

#     with open(os.path.join(path, "inputs-nminors.bin"), "rb") as f:
#         nminors = np.fromfile(f, dtype=np.int32, count=4)
#     nminorlower, nminorklower, nminorupper, nminorkupper = nminors

#     print("Dimensions:", ncol, nlay, nbnd, ngpt, ngas, nflav, neta, npres, ntemp)
#     print("NMinors:", nminorlower, nminorklower, nminorupper, nminorkupper)

#     def read_binary(filename, dtype, shape):
#         """Helper function to read binary file and reshape it."""
#         count = np.prod(shape)
#         data = np.fromfile(os.path.join(path, filename), dtype=dtype, count=count)

#         if shape:
#             return data.reshape(shape, order="C")  # Fortran order reshaping
#         return data

#     data = {
#         "ncol": ncol,
#         "nlay": nlay,
#         "nbnd": nbnd,
#         "ngpt": ngpt,
#         "ngas": ngas,
#         "nflav": nflav,
#         "neta": neta,
#         "npres": npres,
#         "ntemp": ntemp,
#         "nminorlower": nminorlower,
#         "nminorklower": nminorklower,
#         "nminorupper": nminorupper,
#         "nminorkupper": nminorkupper,
#         "idx_h2o": read_binary("inputs-idx_h2o.bin", np.int32, (1,)),
#         "gpoint_flavor": read_binary("inputs-gpoint_flavor.bin", np.int32, (2, ngpt)),
#         "band_lims_gpt": read_binary("inputs-band_lims_gpt.bin", np.int32, (2, nbnd)),
#         "kmajor": read_binary("inputs-kmajor.bin", np.float64, (ntemp, neta, npres+1, ngpt)),
#         "kminor_lower": read_binary("inputs-kminor_lower.bin", np.float64, (ntemp, neta, nminorklower)),
#         "kminor_upper": read_binary("inputs-kminor_upper.bin", np.float64, (ntemp, neta, nminorkupper)),
#         "minor_limits_gpt_lower": read_binary("inputs-minor_limits_gpt_lower.bin", np.int32, (2, nminorlower)),
#         "minor_limits_gpt_upper": read_binary("inputs-minor_limits_gpt_upper.bin", np.int32, (2, nminorupper)),
#         "minor_scales_with_density_lower": read_binary("inputs-minor_scales_with_density_lower.bin", np.bool_, (nminorlower,)),
#         "minor_scales_with_density_upper": read_binary("inputs-minor_scales_with_density_upper.bin", np.bool_, (nminorupper,)),
#         "scale_by_complement_lower": read_binary("inputs-scale_by_complement_lower.bin", np.bool_, (nminorlower,)),
#         "scale_by_complement_upper": read_binary("inputs-scale_by_complement_upper.bin", np.bool_, (nminorupper,)),
#         "idx_minor_lower": read_binary("inputs-idx_minor_lower.bin", np.int32, (nminorlower,)),
#         "idx_minor_upper": read_binary("inputs-idx_minor_upper.bin", np.int32, (nminorupper,)),
#         "idx_minor_scaling_lower": read_binary("inputs-idx_minor_scaling_lower.bin", np.int32, (nminorlower,)),
#         "idx_minor_scaling_upper": read_binary("inputs-idx_minor_scaling_upper.bin", np.int32, (nminorupper,)),
#         "kminor_start_lower": read_binary("inputs-kminor_start_lower.bin", np.int32, (nminorlower,)),
#         "kminor_start_upper": read_binary("inputs-kminor_start_upper.bin", np.int32, (nminorupper,)),
#         "tropo": read_binary("inputs-tropo.bin", np.bool_, (ncol, nlay)),
#         "col_mix": read_binary("inputs-col_mix.bin", np.float64, (2, ncol, nlay, nflav)),
#         "fmajor": read_binary("inputs-fmajor.bin", np.float64, (2, 2, 2, ncol, nlay, nflav)),
#         "fminor": read_binary("inputs-fminor.bin", np.float64, (2, 2, ncol, nlay, nflav)),
#         "play": read_binary("inputs-play.bin", np.float64, (ncol, nlay)),
#         "tlay": read_binary("inputs-tlay.bin", np.float64, (ncol, nlay)),
#         "col_gas": read_binary("inputs-col_gas.bin", np.float64, (ncol, nlay, ngas+1)),
#         "jeta": read_binary("inputs-jeta.bin", np.int32, (2, ncol, nlay, nflav)),
#         "jtemp": read_binary("inputs-jtemp.bin", np.int32, (ncol, nlay)),
#         "jpress": read_binary("inputs-jpress.bin", np.int32, (ncol, nlay)),
#         "tau": read_binary("inputs-tau.bin", np.float64, (ncol, nlay, ngpt))
#     }

#     return data


# def test_rrtmgp_compute_tau_absorption():
#     """
#     Test the Fortran function `test_rrtmgp_compute_tau_absorption` by comparing its outputs
#     against a reference.
#     """
#     output_keys = [
#         "tau"
#     ]

#     # path = "/data/vkm/code/makepath/rte-rrtmgp/build/tests"
#     path = "/data/vkm/code/makepath/rte-rrtmgp/build/examples/all-sky"

#     # Load input data
#     input_data = load_data(path)
#     # # Copy input data to compare outputs
#     output_data = cp.deepcopy(input_data)

#     # Print keys of both dictionaries

#     # Prepare arguments for Fortran function
#     args = list(input_data.values())

#     # Call the Fortran function
#     py.rrtmgp_compute_tau_absorption(*args)

#     # Validate outputs
#     for key in output_keys:
#         ref = output_data[key]
#         computed = input_data[key]

#         diff = np.abs(ref - computed)
#         threshold = 0.1 #7.e-4
#         mask = diff > threshold
#         count = np.sum(mask)
#         if count > 0:
#             indices = np.argwhere(mask)
#             for idx in indices:
#                 ref_val = ref[tuple(idx)]
#                 computed_val = computed[tuple(idx)]
#                 diff_val = diff[tuple(idx)]
#                 # print(f"#{idx}\t{ref_val}\t| {computed_val}\t>> {diff_val}")

#         print(f"{key} shape: {computed.shape}\tDiff count: {count}\tMax diff: {np.max(diff)}")

#         # assert np.allclose(ref, computed, 0.0), f"Mismatch in {key}: {np.max(np.abs(ref - computed))} difference"

# test_rrtmgp_compute_tau_absorption()