"""Definitions of screw axes and space groups for Laue classes."""
from __future__ import absolute_import, division, print_function
from dials.algorithms.symmetry.absences.screw_axes import (
    ScrewAxis21a,
    ScrewAxis21b,
    ScrewAxis21c,
    ScrewAxis41a,
    ScrewAxis41b,
    ScrewAxis41c,
    ScrewAxis42c,
    ScrewAxis61c,
    ScrewAxis62c,
    ScrewAxis63c,
    ScrewAxis31c,
    ScrewAxisObserver,
)


def score_screw_axes(laue_group_info, reflection_table, significance_level=0.95):
    """Get the relevant screw axes and score them. Print pretty table."""
    scores = []
    axes = []
    for axis in laue_group_info["unique_axes"]:
        a = axis()
        a.register_observer(
            event="selected data for scoring", observer=ScrewAxisObserver()
        )
        if "equivalent_axes" in laue_group_info:
            if axis in laue_group_info["equivalent_axes"]:
                for equivalent in laue_group_info["equivalent_axes"][axis]:
                    a.add_equivalent_axis(equivalent())
        scores.append(a.score_axis(reflection_table, significance_level))
        axes.append(a)
    return axes, scores


def score_space_groups(screw_axis_scores, laue_group_info):
    """Score the space groups in a laue group based on axes scores."""
    space_groups = []
    scores = []
    for space_g, conditions in laue_group_info["space_groups"]:
        sg_score = 1.0
        for condition, score in zip(conditions, screw_axis_scores):
            # condition can be True, False or 'pass'
            if condition is True:
                sg_score *= score
            elif condition is False:
                sg_score *= 1.0 - score
        scores.append(sg_score)
        space_groups.append(space_g)
    return space_groups, scores


### For each laue group, define allowed screw axes and how to test these to
### determine the space group; True, False or 'pass' (if not relevant).
### If relevant, list equivalent axes and enantiomorphic pairs.
p2m = {
    "unique_axes": [ScrewAxis21b],
    "space_groups": [("P 2", [False]), ("P 21", [True])],
}

Cmmm = {
    "unique_axes": [ScrewAxis21c],
    "space_groups": [
        ("C 2 2 2", [False]),
        ("C 2 2 21", [True]),  # Assuming we don't allow B 2 21 2, A 21 2 2
    ],
}

Immm = {
    "unique_axes": [ScrewAxis21c, ScrewAxis21b, ScrewAxis21a],
    "space_groups": [
        ("I 2 2 2", [False, False, False]),
        ("I 21 21 21", [True, True, True]),
    ],
}

pmmm = {
    "unique_axes": [ScrewAxis21a, ScrewAxis21b, ScrewAxis21c],
    "space_groups": [
        ("P 2 2 2", [False, False, False]),
        ("P 2 2 21", [False, False, True]),
        ("P 2 21 2", [False, True, False]),
        ("P 21 2 2", [True, False, False]),
        # ^ three equivalent groups, but keep separate to allow standard setting.
        ("P 21 21 2", [True, True, False]),
        ("P 21 2 21", [True, False, True]),
        ("P 2 21 21", [False, True, True]),
        ("P 21 21 21", [True, True, True]),
        # ^ three equivalent groups, but keep separate to allow standard setting.
    ],
}

p6m = {
    "unique_axes": [ScrewAxis61c, ScrewAxis62c, ScrewAxis63c],
    "space_groups": [
        ("P 6", [False, False, False]),
        ("P 61", [True]),
        ("P 62", [False, True]),
        ("P 63", [False, "pass", True]),
    ],
    "enantiomorphic pairs": {"P 61": "P 65", "P 62": "P 64"},
}

p6mmm = {
    "unique_axes": [ScrewAxis61c, ScrewAxis62c, ScrewAxis63c],
    "space_groups": [
        ("P 6 2 2", [False, False, False]),
        ("P 61 2 2", [True]),
        ("P 62 2 2", [False, True]),
        ("P 63 2 2", [False, "pass", True]),
    ],
    "enantiomorphic pairs": {"P 61 2 2": "P 65 2 2", "P 62 2 2": "P 64 2 2"},
}

p3 = {
    "unique_axes": [ScrewAxis31c],
    "space_groups": [("P 3", [False]), ("P 31", [True])],
    "enantiomorphic pairs": {"P 31": "P 32"},
}

p31m = {
    "unique_axes": [ScrewAxis31c],
    "space_groups": [("P 3 1 2", [False]), ("P 31 1 2", [True])],
    "enantiomorphic pairs": {"P 31 1 2": "P 32 1 2"},
}

p3m1 = {
    "unique_axes": [ScrewAxis31c],
    "space_groups": [("P 3 2 1", [False]), ("P 31 2 1", [True])],
    "enantiomorphic pairs": {"P 31 2 1": "P 32 2 1"},
}

p4mmm = {
    "unique_axes": [ScrewAxis41c, ScrewAxis21a, ScrewAxis42c],
    "equivalent_axes": {ScrewAxis21a: [ScrewAxis21b]},
    "space_groups": [
        ("P 4 2 2", [False, False, False]),
        ("P 4 21 2", [False, True, False]),
        ("P 41 2 2", [True, False]),
        ("P 42 2 2", [False, False, True]),
        ("P 41 21 2", [True, True]),
        ("P 42 21 2", [False, True, True]),
    ],
    "enantiomorphic pairs": {"P 41 2 2": "P 43 2 2", "P 41 21 2": "P 43 21 2"},
}

p4m = {
    "unique_axes": [ScrewAxis41c, ScrewAxis42c],
    "space_groups": [
        ("P 4", [False, False]),
        ("P 41", [True]),
        ("P 42", [False, True]),
    ],
    "enantiomorphic pairs": {"P 41": "P 43"},
}

I4m = {
    "unique_axes": [ScrewAxis41c],
    "space_groups": [("I 4", [False]), ("I 41", [True])],
}

I4mmm = {
    "unique_axes": [ScrewAxis41c],
    "space_groups": [("I 4 2 2", [False]), ("I 41 2 2", [True])],
}

pm3 = {
    "unique_axes": [ScrewAxis21a],
    "equivalent_axes": {ScrewAxis21a: [ScrewAxis21b, ScrewAxis21c]},
    "space_groups": [("P 2 3", [False]), ("P 21 3", [True])],
}

Im3 = {
    "unique_axes": [ScrewAxis21a],
    "equivalent_axes": {ScrewAxis21a: [ScrewAxis21b, ScrewAxis21c]},
    "space_groups": [("I 2 3", [False]), ("I 21 3", [True])],
}

pm3m = {
    "unique_axes": [ScrewAxis41a, ScrewAxis42c],
    "equivalent_axes": {ScrewAxis41a: [ScrewAxis41b, ScrewAxis41c]},
    "space_groups": [
        ("P 4 3 2", [False, False]),
        ("P 41 3 2", [True]),
        ("P 42 3 2", [False, True]),
    ],
    "enantiomorphic pairs": {"P 41 3 2": "P 43 3 2"},
}

Im3m = {
    "unique_axes": [ScrewAxis41a],
    "equivalent_axes": {ScrewAxis41a: [ScrewAxis41b, ScrewAxis41c]},
    "space_groups": [("I 4 3 2", [False]), ("I 41 3 2", [True])],
}

Fm3m = {
    "unique_axes": [ScrewAxis41a],
    "equivalent_axes": {ScrewAxis41a: [ScrewAxis41b, ScrewAxis41c]},
    "space_groups": [("F 4 3 2", [False]), ("F 41 3 2", [True])],
}

# These laue groups contain 59 space groups (discarding the equivalent defined in P222)
# - adding P1, C2, F222, R3, R32 and F23 make 65 MX space groups.
laue_groups = {
    "P 1 2/m 1": p2m,
    "P m m m": pmmm,
    "C m m m": Cmmm,
    "I m m m": Immm,
    "P 4/m": p4m,
    "I 4/m": I4m,
    "P 4/m m m": p4mmm,
    "I 4/m m m": I4mmm,
    "P 6/m": p6m,
    "P 6/m m m": p6mmm,
    "P -3": p3,
    "P -3 1 m": p31m,
    "P -3 m 1": p3m1,
    "P m -3": pm3,
    "I m -3": Im3,
    "P m -3 m": pm3m,
    "I m -3 m": Im3m,
    "F m -3 m": Fm3m,
}
