from __future__ import annotations

from scitbx import matrix


def is_approximate_integer_multiple(
    vec_a, vec_b, relative_length_tolerance=0.2, angular_tolerance=5.0
):
    length_a = vec_a.length()
    length_b = vec_b.length()
    # assert length_b >= length_a
    if length_a > length_b:
        vec_a, vec_b = vec_b, vec_a
        length_a, length_b = length_b, length_a
    angle = vec_a.angle(vec_b, deg=True)
    if angle < angular_tolerance or abs(180 - angle) < angular_tolerance:
        n = length_b / length_a
        if abs(round(n) - n) < relative_length_tolerance:
            return True
    return False


class vector_group:
    def __init__(self):
        self.vectors = []
        self.weights = []
        self._mean = None

    def append(self, vector, weight):
        self.vectors.append(vector)
        self.weights.append(weight)
        self._mean = None

    @property
    def mean(self):
        if self._mean is None:
            self._mean = self._compute_mean()
        return self._mean

    def _compute_mean(self):
        sum_x = 0
        sum_y = 0
        sum_z = 0
        for v in self.vectors:
            sum_x += v.elems[0]
            sum_y += v.elems[1]
            sum_z += v.elems[2]
        return matrix.col((sum_x, sum_y, sum_z)) / len(self.vectors)


def group_vectors(
    vectors,
    weights,
    relative_length_tolerance=0.1,
    angular_tolerance=5,
    max_groups=None,
):
    # loop over vectors and sort into groups similar to further down
    # group similar angle and lengths, also catch integer multiples of vectors

    vector_groups = []

    for v, weight in zip(vectors, weights):
        if max_groups is not None and len(vector_groups) == max_groups:
            break
        if not isinstance(v, matrix.col):
            v = matrix.col(v)
        length = v.length()
        matched_group = False
        for group in vector_groups:
            mean_v = group.mean
            mean_v_length = mean_v.length()
            if (
                abs(mean_v_length - length) / max(mean_v_length, length)
                < relative_length_tolerance
            ):
                angle = mean_v.angle(v, deg=True)
                if angle < angular_tolerance:
                    group.append(v, weight)
                    matched_group = True
                    break
                elif abs(180 - angle) < angular_tolerance:
                    group.append(-v, weight)
                    matched_group = True
                    break
        if not matched_group:
            group = vector_group()
            group.append(v, weight)
            vector_groups.append(group)

    return vector_groups
