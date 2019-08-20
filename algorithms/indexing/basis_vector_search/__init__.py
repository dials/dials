from __future__ import absolute_import, division, print_function

from scitbx import matrix


def is_approximate_integer_multiple(
    vec_a, vec_b, relative_tolerance=0.2, angular_tolerance=5.0
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
        if abs(round(n) - n) < relative_tolerance:
            return True
    return False


def find_unique_vectors(vectors, max_vectors):
    unique_vectors = []
    i = 0
    while i < len(vectors) and len(unique_vectors) < max_vectors:
        v = matrix.col(vectors[i])
        is_unique = True
        if i > 0:
            for v_u in unique_vectors:
                if v.length() < v_u.length():
                    if is_approximate_integer_multiple(v, v_u):
                        is_unique = False
                        break
                elif is_approximate_integer_multiple(v_u, v):
                    is_unique = False
                    break
        if is_unique:
            unique_vectors.append(v)
        i += 1
    return unique_vectors
