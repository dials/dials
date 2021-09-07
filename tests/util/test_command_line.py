import random

from libtbx.phil import parse

from dials.util.command_line import OptionParser


def test_guess_input_file_type():
    experiment_files = ["foo.expt", "bar.expt"]
    image_files = ["foo.nxs"] + ["bar_%04d.cbf" % j for j in range(1, 11)]
    reflection_files = ["strong.refl", "indexed.refl"]
    others = ["kitten.jpg", "proposal_final_final_9aug2021.docx"]

    filenames = experiment_files + image_files + reflection_files + others
    random.shuffle(filenames)

    split = OptionParser.guess_input_file_type(filenames)

    assert set(experiment_files) == set(split["experiments"])
    assert set(image_files) == set(split["images"])
    assert set(reflection_files) == set(split["reflections"])
    assert set(others) == set(split["unknown"])


def test_option_parser():
    phil_scope = parse(
        """
    verbose = False
        .type = bool
    egg_price = 197
        .type = int
    """
    )

    op = OptionParser(phil=phil_scope)

    assert op.params().verbose is False
    assert op.params().egg_price == 197
