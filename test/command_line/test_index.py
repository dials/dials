from dxtbx.model.experiment_list import ExperimentList
from dials.array_family import flex
import procrunner


def test_all_expt_ids_have_expts(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.index",
            dials_data("vmxi_thaumatin_grid_index").join("split_07602.expt"),
            dials_data("vmxi_thaumatin_grid_index").join("split_07602.refl"),
            "stills.indexer=sequences",
            "indexing.method=real_space_grid_search",
            "space_group=P4",
            "unit_cell=58,58,150,90,90,90",
            "max_lattices=8",
            "beam.fix=all",
            "detector.fix=all",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("indexed.expt").check(file=1)
    assert tmpdir.join("indexed.refl").check(file=1)

    refl = flex.reflection_table.from_file(tmpdir / "indexed.refl")
    expt = ExperimentList.from_file(tmpdir / "indexed.expt")

    ids = set(refl["id"]) - set([-1])

    assert max(ids) + 1 == len(expt)
