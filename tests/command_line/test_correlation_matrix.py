from __future__ import annotations

import pathlib

import dials.command_line.correlation_matrix as dials_corr_mat


def test_corr_mat(dials_data, run_in_tmp_path):
    mcp = dials_data("vmxi_proteinase_k_sweeps", pathlib=True)
    args = []
    for i in [0, 1, 2, 3]:
        args.append(str(mcp / f"experiments_{i}.expt"))
        args.append(str(mcp / f"reflections_{i}.refl"))
    dials_corr_mat.run(args=args)
    assert pathlib.Path("dials.correlation_matrix.html").is_file()
    assert pathlib.Path("dials.correlation_matrix.log").is_file()
