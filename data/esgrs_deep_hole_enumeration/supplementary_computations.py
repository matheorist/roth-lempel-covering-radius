"""Projective-set equivalence tests for normalized ESGRS extensions.

The implementation is deliberately self-contained and works over prime fields.
It treats nonzero generator-matrix columns as points of PG(r-1, q).  A code
monomial equivalence is therefore recovered as a projectivity between two
unordered projective column sets.

This is a verified small-parameter baseline, not the final recognition
algorithm: its exponent depends on the dimension because it enumerates target
projective frames.
"""

from __future__ import annotations

import argparse
from itertools import combinations, permutations, product
from math import isqrt
from typing import Callable, Iterable, Sequence


Vector = tuple[int, ...]
Matrix = list[list[int]]


def inv_mod(a: int, q: int) -> int:
    a %= q
    if a == 0:
        raise ZeroDivisionError("zero has no inverse")
    return pow(a, -1, q)


def normalize_projective(v: Sequence[int], q: int) -> Vector:
    w = tuple(x % q for x in v)
    for x in w:
        if x:
            z = inv_mod(x, q)
            return tuple((z * y) % q for y in w)
    raise ValueError("the zero vector is not a projective point")


def mat_vec(A: Matrix, v: Sequence[int], q: int) -> Vector:
    return tuple(sum(a * b for a, b in zip(row, v)) % q for row in A)


def mat_mul(A: Matrix, B: Matrix, q: int) -> Matrix:
    rows, inner, cols = len(A), len(B), len(B[0])
    if len(A[0]) != inner:
        raise ValueError("incompatible matrix dimensions")
    return [
        [sum(A[i][t] * B[t][j] for t in range(inner)) % q for j in range(cols)]
        for i in range(rows)
    ]


def mat_inv(A: Matrix, q: int) -> Matrix:
    n = len(A)
    if any(len(row) != n for row in A):
        raise ValueError("matrix must be square")
    aug = [[x % q for x in row] + [int(i == j) for j in range(n)] for i, row in enumerate(A)]
    for col in range(n):
        pivot = next((i for i in range(col, n) if aug[i][col] % q), None)
        if pivot is None:
            raise ValueError("singular matrix")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = inv_mod(aug[col][col], q)
        aug[col] = [(scale * x) % q for x in aug[col]]
        for i in range(n):
            if i != col and aug[i][col] % q:
                factor = aug[i][col] % q
                aug[i] = [(x - factor * y) % q for x, y in zip(aug[i], aug[col])]
    return [row[n:] for row in aug]


def matrix_rank(A: Matrix, q: int) -> int:
    """Rank of a matrix over the prime field F_q."""
    if not A:
        return 0
    B = [[x % q for x in row] for row in A]
    rows, cols = len(B), len(B[0])
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if B[i][col]), None)
        if pivot is None:
            continue
        B[rank], B[pivot] = B[pivot], B[rank]
        scale = inv_mod(B[rank][col], q)
        B[rank] = [(scale * x) % q for x in B[rank]]
        for i in range(rows):
            if i != rank and B[i][col]:
                factor = B[i][col]
                B[i] = [(x - factor * y) % q for x, y in zip(B[i], B[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def row_space_key(rows: Sequence[Sequence[int]], q: int) -> tuple[Vector, ...]:
    """Return the reduced-row-echelon key of a vector subspace."""
    if not rows:
        return ()
    reduced = [[entry % q for entry in row] for row in rows]
    column_count = len(reduced[0])
    if any(len(row) != column_count for row in reduced):
        raise ValueError("rows have inconsistent lengths")
    rank = 0
    for column in range(column_count):
        pivot = next(
            (index for index in range(rank, len(reduced)) if reduced[index][column]),
            None,
        )
        if pivot is None:
            continue
        reduced[rank], reduced[pivot] = reduced[pivot], reduced[rank]
        inverse_pivot = inv_mod(reduced[rank][column], q)
        reduced[rank] = [
            inverse_pivot * entry % q for entry in reduced[rank]
        ]
        for index in range(len(reduced)):
            if index == rank or not reduced[index][column]:
                continue
            factor = reduced[index][column]
            reduced[index] = [
                (entry - factor * pivot_entry) % q
                for entry, pivot_entry in zip(reduced[index], reduced[rank])
            ]
        rank += 1
        if rank == len(reduced):
            break
    return tuple(tuple(row) for row in reduced[:rank])


def columns_to_matrix(columns: Sequence[Sequence[int]]) -> Matrix:
    return [list(row) for row in zip(*columns)]


def is_arc(points: Sequence[Vector], q: int) -> bool:
    """Return whether every r columns of a PG(r-1,q) point set are independent."""
    r = len(points[0])
    return all(
        matrix_rank(columns_to_matrix(tuple(points[i] for i in indices)), q) == r
        for indices in combinations(range(len(points)), r)
    )


def is_frame(points: Sequence[Vector], q: int) -> bool:
    """Return whether r+1 points form a projective frame in PG(r-1,q)."""
    r = len(points[0])
    if len(points) != r + 1:
        return False
    try:
        P_inv = mat_inv(columns_to_matrix(points[:r]), q)
    except ValueError:
        return False
    coordinates = mat_vec(P_inv, points[r], q)
    return all(coordinates)


def frame_projectivity(source: Sequence[Vector], target: Sequence[Vector], q: int) -> Matrix:
    """Construct the unique projectivity carrying one ordered frame to another."""
    r = len(source[0])
    if len(source) != r + 1 or len(target) != r + 1:
        raise ValueError("each projective frame must contain r+1 points")
    P = columns_to_matrix(source[:r])
    Q = columns_to_matrix(target[:r])
    P_inv = mat_inv(P, q)
    Q_inv = mat_inv(Q, q)
    x = mat_vec(P_inv, source[r], q)
    y = mat_vec(Q_inv, target[r], q)
    if not all(x) or not all(y):
        raise ValueError("input is not a projective frame")
    D = [[0] * r for _ in range(r)]
    for i in range(r):
        D[i][i] = y[i] * inv_mod(x[i], q) % q
    return mat_mul(mat_mul(Q, D, q), P_inv, q)


def frame_to_standard_projectivity(frame: Sequence[Vector], q: int) -> Matrix:
    """Map an ordered projective frame to e_0,...,e_{r-1},(1,...,1)."""
    r = len(frame[0])
    if len(frame) != r + 1:
        raise ValueError("a projective frame in PG(r-1,q) has r+1 points")
    P_inv = mat_inv(columns_to_matrix(frame[:r]), q)
    x = mat_vec(P_inv, frame[r], q)
    if not all(x):
        raise ValueError("input is not a projective frame")
    D = [[0] * r for _ in range(r)]
    for i in range(r):
        D[i][i] = inv_mod(x[i], q)
    return mat_mul(D, P_inv, q)


def first_frame(points: Sequence[Vector], q: int) -> tuple[Vector, ...]:
    r = len(points[0])
    for indices in permutations(range(len(points)), r + 1):
        candidate = tuple(points[i] for i in indices)
        if is_frame(candidate, q):
            return candidate
    raise ValueError("the point set contains no projective frame")


def projective_equivalence(
    source: Iterable[Sequence[int]], target: Iterable[Sequence[int]], q: int
) -> Matrix | None:
    """Return an explicit projectivity between two unordered point sets, if one exists."""
    src = tuple(normalize_projective(v, q) for v in source)
    dst = tuple(normalize_projective(v, q) for v in target)
    if len(src) != len(dst) or len(set(src)) != len(src) or len(set(dst)) != len(dst):
        raise ValueError("point sets must have the same size and contain no repetitions")
    source_frame = first_frame(src, q)
    target_set = set(dst)
    r = len(src[0])
    for indices in permutations(range(len(dst)), r + 1):
        target_frame = tuple(dst[i] for i in indices)
        if not is_frame(target_frame, q):
            continue
        T = frame_projectivity(source_frame, target_frame, q)
        image = {normalize_projective(mat_vec(T, p, q), q) for p in src}
        if image == target_set:
            return T
    return None


def monomial_map_from_projectivity(
    source: Sequence[Vector], target: Sequence[Vector], projectivity: Matrix, q: int
) -> dict[str, tuple[int, ...]]:
    """Recover the right monomial map induced by an ambient projectivity.

    If ``T`` is ``projectivity`` and the result is ``(sigma, scales)``, then
    target coordinate ``j`` receives source coordinate ``sigma[j]`` multiplied
    by ``scales[j]``.  Consequently ``T * G_source * M = G_target`` exactly,
    not merely projectively.
    """
    if len(source) != len(target):
        raise ValueError("source and target must have the same number of columns")
    target_by_point = {
        normalize_projective(point, q): index for index, point in enumerate(target)
    }
    if len(target_by_point) != len(target):
        raise ValueError("target projective columns must be distinct")
    source_for_target = [-1] * len(target)
    target_from_source = [-1] * len(source)
    right_scales = [0] * len(target)
    for source_index, source_point in enumerate(source):
        transformed = mat_vec(projectivity, source_point, q)
        key = normalize_projective(transformed, q)
        if key not in target_by_point:
            raise ValueError("projectivity does not map the source set to the target set")
        target_index = target_by_point[key]
        if source_for_target[target_index] != -1:
            raise ValueError("projectivity does not induce a coordinate permutation")
        target_point = tuple(value % q for value in target[target_index])
        pivot = next(index for index, value in enumerate(target_point) if value)
        projective_scale = transformed[pivot] * inv_mod(target_point[pivot], q) % q
        right_scale = inv_mod(projective_scale, q)
        assert tuple(value * right_scale % q for value in transformed) == target_point
        source_for_target[target_index] = source_index
        target_from_source[source_index] = target_index
        right_scales[target_index] = right_scale
    assert -1 not in source_for_target and -1 not in target_from_source
    return {
        "source_for_target": tuple(source_for_target),
        "target_from_source": tuple(target_from_source),
        "right_scales": tuple(right_scales),
    }


def apply_right_monomial(
    word: Sequence[int], monomial_map: dict[str, tuple[int, ...]], q: int, *, inverse: bool = False
) -> Vector:
    """Apply the recovered right monomial map, or its inverse, to a row word."""
    source_for_target = monomial_map["source_for_target"]
    right_scales = monomial_map["right_scales"]
    if len(word) != len(source_for_target):
        raise ValueError("word length does not match the monomial map")
    output = [0] * len(word)
    if not inverse:
        for target_index, source_index in enumerate(source_for_target):
            output[target_index] = word[source_index] * right_scales[target_index] % q
    else:
        for target_index, source_index in enumerate(source_for_target):
            output[source_index] = (
                word[target_index] * inv_mod(right_scales[target_index], q) % q
            )
    return tuple(output)


def transport_unique_decoder(
    received: Sequence[int],
    monomial_map: dict[str, tuple[int, ...]],
    target_decoder: Callable[[Sequence[int]], Sequence[int]],
    q: int,
) -> Vector:
    """Decode through a recovered monomial equivalence and transport back."""
    target_received = apply_right_monomial(received, monomial_map, q)
    target_codeword = target_decoder(target_received)
    return apply_right_monomial(target_codeword, monomial_map, q, inverse=True)


def transport_list_decoder(
    received: Sequence[int],
    monomial_map: dict[str, tuple[int, ...]],
    target_decoder: Callable[[Sequence[int]], Iterable[Sequence[int]]],
    q: int,
) -> tuple[Vector, ...]:
    """List decode through a recovered monomial equivalence and transport back."""
    target_received = apply_right_monomial(received, monomial_map, q)
    target_codewords = target_decoder(target_received)
    return tuple(
        apply_right_monomial(codeword, monomial_map, q, inverse=True)
        for codeword in target_codewords
    )


def encode_from_columns(message: Sequence[int], columns: Sequence[Vector], q: int) -> Vector:
    """Encode a row message using a generator matrix given by its columns."""
    if not columns or len(message) != len(columns[0]):
        raise ValueError("message dimension does not match the generator columns")
    return tuple(
        sum(coefficient * entry for coefficient, entry in zip(message, column)) % q
        for column in columns
    )


def right_nullspace_basis(A: Matrix, q: int) -> tuple[Vector, ...]:
    """Return a row basis for the right nullspace of ``A`` over F_q."""
    if not A or not A[0]:
        raise ValueError("the matrix must be nonempty")
    width = len(A[0])
    if any(len(row) != width for row in A):
        raise ValueError("matrix rows must have a common length")
    B = [[entry % q for entry in row] for row in A]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(width):
        pivot = next((row for row in range(pivot_row, len(B)) if B[row][column]), None)
        if pivot is None:
            continue
        B[pivot_row], B[pivot] = B[pivot], B[pivot_row]
        scale = inv_mod(B[pivot_row][column], q)
        B[pivot_row] = [(scale * entry) % q for entry in B[pivot_row]]
        for row in range(len(B)):
            if row != pivot_row and B[row][column]:
                factor = B[row][column]
                B[row] = [
                    (entry - factor * pivot_entry) % q
                    for entry, pivot_entry in zip(B[row], B[pivot_row])
                ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(B):
            break

    free_columns = [column for column in range(width) if column not in pivot_columns]
    basis = []
    for free_column in free_columns:
        vector = [0] * width
        vector[free_column] = 1
        for row, pivot_column in enumerate(pivot_columns):
            vector[pivot_column] = -B[row][free_column] % q
        basis.append(tuple(vector))
    return tuple(basis)


def parity_check_from_columns(columns: Sequence[Vector], q: int) -> tuple[Vector, ...]:
    """Construct a full-row-rank parity-check matrix from generator columns."""
    if not columns:
        raise ValueError("the generator matrix must have at least one column")
    generator = columns_to_matrix(columns)
    dimension = len(generator)
    if matrix_rank(generator, q) != dimension:
        raise ValueError("the generator matrix does not have full row rank")
    parity_check = right_nullspace_basis(generator, q)
    if len(parity_check) != len(columns) - dimension:
        raise AssertionError("unexpected parity-check dimension")
    assert all(
        sum(generator[row][column] * parity_check[check][column] for column in range(len(columns)))
        % q
        == 0
        for row in range(dimension)
        for check in range(len(parity_check))
    )
    return parity_check


def word_syndrome(received: Sequence[int], parity_check: Sequence[Vector], q: int) -> Vector:
    """Compute H received^T when the rows of H are given by ``parity_check``."""
    if not parity_check or len(received) != len(parity_check[0]):
        raise ValueError("received word and parity-check matrix have incompatible lengths")
    return tuple(
        sum(entry * symbol for entry, symbol in zip(row, received)) % q
        for row in parity_check
    )


def prepare_single_error_syndrome_decoder(
    columns: Sequence[Vector], q: int
) -> Callable[[Sequence[int]], Vector]:
    """Precompute a linear-time decoder for a code of minimum distance at least three.

    A nonzero one-error syndrome is a scalar multiple of exactly one parity-check
    column.  The projective lookup table therefore locates the error coordinate,
    after which one division recovers its value.  For an MDS [n,n-3,4] code this
    reaches the optimal unique-decoding radius.
    """
    parity_check = parity_check_from_columns(columns, q)
    redundancy = len(parity_check)
    check_columns = tuple(
        tuple(parity_check[row][column] for row in range(redundancy))
        for column in range(len(columns))
    )
    locator: dict[Vector, tuple[int, Vector]] = {}
    for coordinate, check_column in enumerate(check_columns):
        if not any(check_column):
            raise ValueError("a zero parity-check column prevents single-error correction")
        key = normalize_projective(check_column, q)
        if key in locator:
            raise ValueError("proportional parity-check columns prevent unique error location")
        locator[key] = (coordinate, check_column)

    def decode(received: Sequence[int]) -> Vector:
        if len(received) != len(columns):
            raise ValueError("received word has the wrong length")
        normalized_received = tuple(symbol % q for symbol in received)
        syndrome = word_syndrome(normalized_received, parity_check, q)
        if not any(syndrome):
            return normalized_received
        location = locator.get(normalize_projective(syndrome, q))
        if location is None:
            raise ValueError("the syndrome is not consistent with at most one error")
        coordinate, check_column = location
        pivot = next(index for index, entry in enumerate(check_column) if entry)
        error_value = syndrome[pivot] * inv_mod(check_column[pivot], q) % q
        corrected = list(normalized_received)
        corrected[coordinate] = (corrected[coordinate] - error_value) % q
        if any(word_syndrome(corrected, parity_check, q)):
            raise AssertionError("single-error correction did not clear the syndrome")
        return tuple(corrected)

    return decode


def quotient_projection_with_kernel(vector: Sequence[int], q: int) -> Matrix:
    """Return a rank-(r-1) matrix whose kernel is the span of an r-vector."""
    if len(vector) < 2 or not any(entry % q for entry in vector):
        raise ValueError("a nonzero vector of length at least two is required")
    v = tuple(entry % q for entry in vector)
    pivot = next(index for index, entry in enumerate(v) if entry)
    rows: Matrix = []
    for index in range(len(v)):
        if index == pivot:
            continue
        row = [0] * len(v)
        row[pivot] = -v[index] % q
        row[index] = v[pivot]
        rows.append(row)
    assert matrix_rank(rows, q) == len(v) - 1
    assert mat_vec(rows, v, q) == tuple(0 for _ in range(len(v) - 1))
    return rows


def solve_independent_column_combination(
    vectors: Sequence[Sequence[int]], target: Sequence[int], q: int
) -> tuple[int, ...]:
    """Solve sum(c_i*v_i)=target when the input columns are independent."""
    if not vectors:
        raise ValueError("at least one column is required")
    ambient_dimension = len(target)
    column_count = len(vectors)
    if column_count > ambient_dimension or any(
        len(vector) != ambient_dimension for vector in vectors
    ):
        raise ValueError("column dimensions are incompatible with the target")
    for row_indices in combinations(range(ambient_dimension), column_count):
        minor = [
            [vectors[column][row] % q for column in range(column_count)]
            for row in row_indices
        ]
        try:
            inverse_minor = mat_inv(minor, q)
        except ValueError:
            continue
        coefficients = mat_vec(inverse_minor, tuple(target[row] for row in row_indices), q)
        reconstructed = tuple(
            sum(coefficient * vector[row] for coefficient, vector in zip(coefficients, vectors))
            % q
            for row in range(ambient_dimension)
        )
        if reconstructed != tuple(entry % q for entry in target):
            raise ValueError("the target is not in the span of the columns")
        return coefficients
    raise ValueError("the input columns are linearly dependent")


def solve_two_column_combination(
    left: Sequence[int], right: Sequence[int], target: Sequence[int], q: int
) -> tuple[int, int]:
    """Solve a*left+b*right=target for two independent columns."""
    coefficients = solve_independent_column_combination((left, right), target, q)
    return coefficients[0], coefficients[1]


def projective_line_key(left: Sequence[int], right: Sequence[int], q: int) -> Vector:
    """Return normalized Pluecker coordinates of a projective line."""
    if len(left) < 3 or len(left) != len(right):
        raise ValueError("projective points require equal dimension at least three")
    line = tuple(
        (left[i] * right[j] - left[j] * right[i]) % q
        for i in range(len(left))
        for j in range(i + 1, len(left))
    )
    return normalize_projective(line, q)


def projective_plane_key(
    first: Sequence[int], second: Sequence[int], third: Sequence[int], q: int
) -> Vector:
    """Return the normalized dual coordinates of a plane in PG(3,q)."""
    rows = tuple(tuple(entry % q for entry in vector) for vector in (first, second, third))
    if any(len(row) != 4 for row in rows):
        raise ValueError("projective three-space points require four coordinates")
    if matrix_rank([list(row) for row in rows], q) != 3:
        raise ValueError("three noncollinear projective points are required")

    def determinant_three(matrix: Sequence[Sequence[int]]) -> int:
        return (
            matrix[0][0]
            * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
            - matrix[0][1]
            * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
            + matrix[0][2]
            * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
        ) % q

    plane = []
    for omitted_column in range(4):
        minor = [
            [row[column] for column in range(4) if column != omitted_column]
            for row in rows
        ]
        cofactor = determinant_three(minor)
        plane.append(cofactor if omitted_column % 2 == 0 else -cofactor % q)
    key = normalize_projective(plane, q)
    assert all(sum(x * y for x, y in zip(row, key)) % q == 0 for row in rows)
    return key


def prepare_radius_two_syndrome_list_decoder(
    columns: Sequence[Vector], q: int
) -> Callable[[Sequence[int]], tuple[Vector, ...]]:
    """Precompute a linear-time radius-two list decoder for an MDS [n,n-3,4] code.

    For nonzero syndrome ``s``, quotienting the three-dimensional syndrome
    space by ``<s>`` turns every secant through ``[s]`` into a collision in
    PG(1,q).  Hashing these projected parity-check columns finds all weight-two
    error patterns in linear time.  The arc property bounds every bucket by two.
    """
    parity_check = parity_check_from_columns(columns, q)
    if len(parity_check) != 3:
        raise ValueError("the radius-two quotient decoder requires redundancy three")
    check_columns = tuple(
        tuple(parity_check[row][column] for row in range(3))
        for column in range(len(columns))
    )
    projective_columns = [normalize_projective(column, q) for column in check_columns]
    if len(set(projective_columns)) != len(projective_columns):
        raise ValueError("the parity-check columns do not define a projective arc")

    def decode(received: Sequence[int]) -> tuple[Vector, ...]:
        if len(received) != len(columns):
            raise ValueError("received word has the wrong length")
        normalized_received = tuple(symbol % q for symbol in received)
        syndrome = word_syndrome(normalized_received, parity_check, q)
        if not any(syndrome):
            return (normalized_received,)

        quotient = quotient_projection_with_kernel(syndrome, q)
        buckets: dict[Vector, list[int]] = {}
        one_error_coordinate: int | None = None
        for coordinate, check_column in enumerate(check_columns):
            projection = mat_vec(quotient, check_column, q)
            if not any(projection):
                if one_error_coordinate is not None:
                    raise ValueError("multiple check columns are proportional to the syndrome")
                one_error_coordinate = coordinate
                continue
            key = normalize_projective(projection, q)
            bucket = buckets.setdefault(key, [])
            bucket.append(coordinate)
            if len(bucket) > 2:
                raise ValueError("three check columns lie on a line through the syndrome")

        candidates: set[Vector] = set()
        if one_error_coordinate is not None:
            check_column = check_columns[one_error_coordinate]
            pivot = next(index for index, entry in enumerate(check_column) if entry)
            error_value = syndrome[pivot] * inv_mod(check_column[pivot], q) % q
            corrected = list(normalized_received)
            corrected[one_error_coordinate] = (
                corrected[one_error_coordinate] - error_value
            ) % q
            candidates.add(tuple(corrected))

        for bucket in buckets.values():
            if len(bucket) != 2:
                continue
            left_coordinate, right_coordinate = bucket
            left_error, right_error = solve_two_column_combination(
                check_columns[left_coordinate],
                check_columns[right_coordinate],
                syndrome,
                q,
            )
            if not left_error or not right_error:
                raise AssertionError("a secant collision produced a zero error value")
            corrected = list(normalized_received)
            corrected[left_coordinate] = (
                corrected[left_coordinate] - left_error
            ) % q
            corrected[right_coordinate] = (
                corrected[right_coordinate] - right_error
            ) % q
            candidates.add(tuple(corrected))

        assert all(not any(word_syndrome(candidate, parity_check, q)) for candidate in candidates)
        assert all(
            sum(x != y for x, y in zip(candidate, normalized_received)) <= 2
            for candidate in candidates
        )
        return tuple(sorted(candidates))

    return decode


def prepare_radius_three_syndrome_list_decoder(
    columns: Sequence[Vector],
    q: int,
    *,
    parity_check: Sequence[Sequence[int]] | None = None,
) -> Callable[[Sequence[int]], tuple[Vector, ...]]:
    """Precompute a radius-three decoder for an MDS [n,n-R,R+1] code, R>=4.

    Quotienting the R-dimensional syndrome space by ``<s>`` maps exact
    weight-three supports to collinear triples in PG(R-2,q).  Pluecker-line
    hashing finds all such triples in expected O(n^2) time for fixed R, one
    exponent below direct enumeration of coordinate triples.
    """
    if parity_check is None:
        check_rows = parity_check_from_columns(columns, q)
    else:
        check_rows = tuple(tuple(entry % q for entry in row) for row in parity_check)
        if not check_rows or any(len(row) != len(columns) for row in check_rows):
            raise ValueError("parity-check rows have the wrong length")
        generator = columns_to_matrix(columns)
        if any(
            sum(generator[row][coordinate] * check_rows[check][coordinate]
                for coordinate in range(len(columns))) % q
            for row in range(len(generator))
            for check in range(len(check_rows))
        ):
            raise ValueError("the supplied parity-check matrix is not orthogonal")
    redundancy = len(check_rows)
    if redundancy < 4 or matrix_rank([list(row) for row in check_rows], q) != redundancy:
        raise ValueError("the radius-three quotient decoder requires redundancy at least four")
    check_columns = tuple(
        tuple(check_rows[row][column] for row in range(redundancy))
        for column in range(len(columns))
    )
    projective_columns = [normalize_projective(column, q) for column in check_columns]
    if len(set(projective_columns)) != len(projective_columns):
        raise ValueError("the parity-check columns do not define a projective arc")

    def decode(received: Sequence[int]) -> tuple[Vector, ...]:
        if len(received) != len(columns):
            raise ValueError("received word has the wrong length")
        normalized_received = tuple(symbol % q for symbol in received)
        syndrome = word_syndrome(normalized_received, check_rows, q)
        if not any(syndrome):
            return (normalized_received,)

        quotient = quotient_projection_with_kernel(syndrome, q)
        projection_buckets: dict[Vector, list[int]] = {}
        one_error_coordinate: int | None = None
        for coordinate, check_column in enumerate(check_columns):
            projection = mat_vec(quotient, check_column, q)
            if not any(projection):
                if one_error_coordinate is not None:
                    raise ValueError("multiple check columns are proportional to the syndrome")
                one_error_coordinate = coordinate
                continue
            key = normalize_projective(projection, q)
            bucket = projection_buckets.setdefault(key, [])
            bucket.append(coordinate)
            if len(bucket) > 2:
                raise ValueError("three check columns span a plane containing the syndrome")

        candidates: set[Vector] = set()
        if one_error_coordinate is not None:
            check_column = check_columns[one_error_coordinate]
            pivot = next(index for index, entry in enumerate(check_column) if entry)
            error_value = syndrome[pivot] * inv_mod(check_column[pivot], q) % q
            corrected = list(normalized_received)
            corrected[one_error_coordinate] = (
                corrected[one_error_coordinate] - error_value
            ) % q
            candidates.add(tuple(corrected))

        unique_projection_coordinates: list[int] = []
        projections_by_coordinate: dict[int, Vector] = {}
        for key, bucket in projection_buckets.items():
            if len(bucket) == 1:
                coordinate = bucket[0]
                unique_projection_coordinates.append(coordinate)
                projections_by_coordinate[coordinate] = key
                continue
            left_coordinate, right_coordinate = bucket
            errors = solve_independent_column_combination(
                (check_columns[left_coordinate], check_columns[right_coordinate]),
                syndrome,
                q,
            )
            if not all(errors):
                raise AssertionError("a secant collision produced a zero error value")
            corrected = list(normalized_received)
            for coordinate, error_value in zip(bucket, errors):
                corrected[coordinate] = (corrected[coordinate] - error_value) % q
            candidates.add(tuple(corrected))

        line_buckets: dict[Vector, set[int]] = {}
        for left_coordinate, right_coordinate in combinations(
            unique_projection_coordinates, 2
        ):
            line = projective_line_key(
                projections_by_coordinate[left_coordinate],
                projections_by_coordinate[right_coordinate],
                q,
            )
            coordinates_on_line = line_buckets.setdefault(line, set())
            coordinates_on_line.update((left_coordinate, right_coordinate))
            if len(coordinates_on_line) > 3:
                raise ValueError("four projected check columns are collinear")

        for coordinates_on_line in line_buckets.values():
            if len(coordinates_on_line) != 3:
                continue
            support = tuple(sorted(coordinates_on_line))
            errors = solve_independent_column_combination(
                tuple(check_columns[coordinate] for coordinate in support),
                syndrome,
                q,
            )
            if not all(errors):
                raise AssertionError("a trisecant collision produced a zero error value")
            corrected = list(normalized_received)
            for coordinate, error_value in zip(support, errors):
                corrected[coordinate] = (corrected[coordinate] - error_value) % q
            candidates.add(tuple(corrected))

        assert all(not any(word_syndrome(candidate, check_rows, q)) for candidate in candidates)
        assert all(
            sum(x != y for x, y in zip(candidate, normalized_received)) <= 3
            for candidate in candidates
        )
        return tuple(sorted(candidates))

    return decode


def prepare_radius_four_syndrome_list_decoder(
    columns: Sequence[Vector],
    q: int,
    *,
    parity_check: Sequence[Sequence[int]] | None = None,
) -> Callable[[Sequence[int]], tuple[Vector, ...]]:
    """Precompute a cubic-time radius-four decoder for an MDS [n,n-5,6] code.

    Quotienting syndrome space by ``<s>`` maps exact weight-four supports to
    four-point projective circuits in PG(3,q).  Hashing the plane through each
    noncollinear projected triple finds all such supports in expected O(n^3)
    time.
    """
    if parity_check is None:
        check_rows = parity_check_from_columns(columns, q)
    else:
        check_rows = tuple(tuple(entry % q for entry in row) for row in parity_check)
        if not check_rows or any(len(row) != len(columns) for row in check_rows):
            raise ValueError("parity-check rows have the wrong length")
        generator = columns_to_matrix(columns)
        if any(
            sum(
                generator[row][coordinate] * check_rows[check][coordinate]
                for coordinate in range(len(columns))
            )
            % q
            for row in range(len(generator))
            for check in range(len(check_rows))
        ):
            raise ValueError("the supplied parity-check matrix is not orthogonal")
    if len(check_rows) != 5 or matrix_rank([list(row) for row in check_rows], q) != 5:
        raise ValueError("the radius-four quotient decoder requires redundancy five")
    check_columns = tuple(
        tuple(check_rows[row][column] for row in range(5))
        for column in range(len(columns))
    )
    projective_columns = [normalize_projective(column, q) for column in check_columns]
    if len(set(projective_columns)) != len(projective_columns):
        raise ValueError("the parity-check columns do not define a projective arc")

    def decode(received: Sequence[int]) -> tuple[Vector, ...]:
        if len(received) != len(columns):
            raise ValueError("received word has the wrong length")
        normalized_received = tuple(symbol % q for symbol in received)
        syndrome = word_syndrome(normalized_received, check_rows, q)
        if not any(syndrome):
            return (normalized_received,)

        quotient = quotient_projection_with_kernel(syndrome, q)
        projection_buckets: dict[Vector, list[int]] = {}
        one_error_coordinate: int | None = None
        for coordinate, check_column in enumerate(check_columns):
            projection = mat_vec(quotient, check_column, q)
            if not any(projection):
                if one_error_coordinate is not None:
                    raise ValueError("multiple check columns are proportional to the syndrome")
                one_error_coordinate = coordinate
                continue
            key = normalize_projective(projection, q)
            bucket = projection_buckets.setdefault(key, [])
            bucket.append(coordinate)
            if len(bucket) > 2:
                raise ValueError("three check columns span a plane containing the syndrome")

        candidates: set[Vector] = set()

        def record_candidate(support: Sequence[int]) -> None:
            errors = solve_independent_column_combination(
                tuple(check_columns[coordinate] for coordinate in support),
                syndrome,
                q,
            )
            if not all(errors):
                return
            corrected = list(normalized_received)
            for coordinate, error_value in zip(support, errors):
                corrected[coordinate] = (corrected[coordinate] - error_value) % q
            candidates.add(tuple(corrected))

        if one_error_coordinate is not None:
            record_candidate((one_error_coordinate,))

        unique_projection_coordinates: list[int] = []
        projections_by_coordinate: dict[int, Vector] = {}
        for key, bucket in projection_buckets.items():
            if len(bucket) == 1:
                coordinate = bucket[0]
                unique_projection_coordinates.append(coordinate)
                projections_by_coordinate[coordinate] = key
            else:
                record_candidate(tuple(bucket))

        line_buckets: dict[Vector, set[int]] = {}
        for left_coordinate, right_coordinate in combinations(
            unique_projection_coordinates, 2
        ):
            line = projective_line_key(
                projections_by_coordinate[left_coordinate],
                projections_by_coordinate[right_coordinate],
                q,
            )
            coordinates_on_line = line_buckets.setdefault(line, set())
            coordinates_on_line.update((left_coordinate, right_coordinate))
            if len(coordinates_on_line) > 3:
                raise ValueError("four projected check columns are collinear")
        for coordinates_on_line in line_buckets.values():
            if len(coordinates_on_line) == 3:
                record_candidate(tuple(sorted(coordinates_on_line)))

        plane_buckets: dict[Vector, set[int]] = {}
        for support in combinations(unique_projection_coordinates, 3):
            projected = tuple(projections_by_coordinate[coordinate] for coordinate in support)
            if matrix_rank(columns_to_matrix(projected), q) != 3:
                continue
            plane = projective_plane_key(*projected, q)
            coordinates_on_plane = plane_buckets.setdefault(plane, set())
            coordinates_on_plane.update(support)
            if len(coordinates_on_plane) > 4:
                raise ValueError("five projected check columns are coplanar")
        for coordinates_on_plane in plane_buckets.values():
            if len(coordinates_on_plane) == 4:
                record_candidate(tuple(sorted(coordinates_on_plane)))

        assert all(not any(word_syndrome(candidate, check_rows, q)) for candidate in candidates)
        assert all(
            sum(x != y for x, y in zip(candidate, normalized_received)) <= 4
            for candidate in candidates
        )
        return tuple(sorted(candidates))

    return decode


def verify_all_radius_two_error_patterns(
    columns: Sequence[Vector],
    codeword: Sequence[int],
    list_decoder: Callable[[Sequence[int]], Iterable[Sequence[int]]],
    q: int,
) -> tuple[int, int]:
    """Regression-test a list decoder on every weight-at-most-two error pattern."""
    parity_check = parity_check_from_columns(columns, q)
    test_count = 0
    maximum_list_size = 0
    error_patterns = [()]
    error_patterns.extend(
        ((coordinate, error_value),)
        for coordinate in range(len(codeword))
        for error_value in range(1, q)
    )
    error_patterns.extend(
        ((left, left_value), (right, right_value))
        for left, right in combinations(range(len(codeword)), 2)
        for left_value in range(1, q)
        for right_value in range(1, q)
    )
    for error_pattern in error_patterns:
        received = list(codeword)
        for coordinate, error_value in error_pattern:
            received[coordinate] = (received[coordinate] + error_value) % q
        decoded_list = tuple(list_decoder(received))
        assert tuple(codeword) in decoded_list
        assert len(decoded_list) <= len(codeword) // 2
        assert all(
            not any(word_syndrome(candidate, parity_check, q))
            and sum(x != y for x, y in zip(candidate, received)) <= 2
            for candidate in decoded_list
        )
        test_count += 1
        maximum_list_size = max(maximum_list_size, len(decoded_list))
    return test_count, maximum_list_size


def verify_all_radius_two_syndrome_cosets(
    columns: Sequence[Vector],
    list_decoder: Callable[[Sequence[int]], Iterable[Sequence[int]]],
    q: int,
) -> tuple[int, int]:
    """Compare a radius-two decoder with all weight-at-most-two errors in every coset."""
    parity_check = parity_check_from_columns(columns, q)
    if len(parity_check) != 3:
        raise ValueError("the coset regression requires redundancy three")
    check_columns = tuple(
        tuple(parity_check[row][column] for row in range(3))
        for column in range(len(columns))
    )
    basis_indices = next(
        indices
        for indices in combinations(range(len(columns)), 3)
        if matrix_rank(columns_to_matrix(tuple(check_columns[i] for i in indices)), q) == 3
    )
    basis_inverse = mat_inv(
        columns_to_matrix(tuple(check_columns[i] for i in basis_indices)), q
    )

    errors_by_syndrome: dict[Vector, set[Vector]] = {}

    def record_error(entries: Sequence[tuple[int, int]]) -> None:
        error = [0] * len(columns)
        for coordinate, value in entries:
            error[coordinate] = value
        syndrome = word_syndrome(error, parity_check, q)
        errors_by_syndrome.setdefault(syndrome, set()).add(tuple(error))

    record_error(())
    for coordinate in range(len(columns)):
        for value in range(1, q):
            record_error(((coordinate, value),))
    for left, right in combinations(range(len(columns)), 2):
        for left_value in range(1, q):
            for right_value in range(1, q):
                record_error(((left, left_value), (right, right_value)))

    maximum_list_size = 0
    coset_count = 0
    for syndrome in product(range(q), repeat=3):
        received = [0] * len(columns)
        basis_values = mat_vec(basis_inverse, syndrome, q)
        for coordinate, value in zip(basis_indices, basis_values):
            received[coordinate] = value
        assert word_syndrome(received, parity_check, q) == syndrome
        decoded_list = tuple(list_decoder(received))
        recovered_errors = {
            tuple((symbol - decoded_symbol) % q for symbol, decoded_symbol in zip(received, candidate))
            for candidate in decoded_list
        }
        assert recovered_errors == errors_by_syndrome.get(syndrome, set())
        maximum_list_size = max(maximum_list_size, len(decoded_list))
        coset_count += 1
    return coset_count, maximum_list_size


def verify_all_radius_three_syndrome_cosets(
    columns: Sequence[Vector],
    list_decoder: Callable[[Sequence[int]], Iterable[Sequence[int]]],
    q: int,
    *,
    parity_check: Sequence[Sequence[int]] | None = None,
) -> tuple[int, int, dict[int, int], int]:
    """Verify exact radius-three output lists on every redundancy-four coset."""
    check_rows = (
        parity_check_from_columns(columns, q)
        if parity_check is None
        else tuple(tuple(entry % q for entry in row) for row in parity_check)
    )
    if len(check_rows) != 4:
        raise ValueError("the coset regression requires redundancy four")
    check_columns = tuple(
        tuple(check_rows[row][column] for row in range(4))
        for column in range(len(columns))
    )
    basis_indices = next(
        indices
        for indices in combinations(range(len(columns)), 4)
        if matrix_rank(columns_to_matrix(tuple(check_columns[i] for i in indices)), q) == 4
    )
    basis_inverse = mat_inv(
        columns_to_matrix(tuple(check_columns[i] for i in basis_indices)), q
    )

    errors_by_syndrome: dict[Vector, set[Vector]] = {}

    def record_error(entries: Sequence[tuple[int, int]]) -> None:
        error = [0] * len(columns)
        for coordinate, value in entries:
            error[coordinate] = value
        syndrome = word_syndrome(error, check_rows, q)
        errors_by_syndrome.setdefault(syndrome, set()).add(tuple(error))

    record_error(())
    for coordinate in range(len(columns)):
        for value in range(1, q):
            record_error(((coordinate, value),))
    for left, right in combinations(range(len(columns)), 2):
        for left_value in range(1, q):
            for right_value in range(1, q):
                record_error(((left, left_value), (right, right_value)))
    for support in combinations(range(len(columns)), 3):
        for values in product(range(1, q), repeat=3):
            record_error(tuple(zip(support, values)))

    list_size_histogram: dict[int, int] = {}
    maximum_list_size = 0
    coset_count = 0
    for syndrome in product(range(q), repeat=4):
        received = [0] * len(columns)
        basis_values = mat_vec(basis_inverse, syndrome, q)
        for coordinate, value in zip(basis_indices, basis_values):
            received[coordinate] = value
        assert word_syndrome(received, check_rows, q) == syndrome
        decoded_list = tuple(list_decoder(received))
        recovered_errors = {
            tuple((symbol - decoded_symbol) % q for symbol, decoded_symbol in zip(received, candidate))
            for candidate in decoded_list
        }
        assert recovered_errors == errors_by_syndrome.get(syndrome, set())
        list_size = len(decoded_list)
        list_size_histogram[list_size] = list_size_histogram.get(list_size, 0) + 1
        maximum_list_size = max(maximum_list_size, list_size)
        coset_count += 1
    error_count = sum(len(errors) for errors in errors_by_syndrome.values())
    return coset_count, maximum_list_size, list_size_histogram, error_count


def verify_all_radius_four_syndrome_cosets(
    columns: Sequence[Vector],
    list_decoder: Callable[[Sequence[int]], Iterable[Sequence[int]]],
    q: int,
    *,
    parity_check: Sequence[Sequence[int]] | None = None,
) -> tuple[int, int, dict[int, int], int]:
    """Verify exact radius-four lists on every redundancy-five coset."""
    check_rows = (
        parity_check_from_columns(columns, q)
        if parity_check is None
        else tuple(tuple(entry % q for entry in row) for row in parity_check)
    )
    if len(check_rows) != 5:
        raise ValueError("the coset regression requires redundancy five")
    check_columns = tuple(
        tuple(check_rows[row][column] for row in range(5))
        for column in range(len(columns))
    )
    basis_indices = next(
        indices
        for indices in combinations(range(len(columns)), 5)
        if matrix_rank(columns_to_matrix(tuple(check_columns[i] for i in indices)), q) == 5
    )
    basis_inverse = mat_inv(
        columns_to_matrix(tuple(check_columns[i] for i in basis_indices)), q
    )

    errors_by_syndrome: dict[Vector, set[Vector]] = {}

    def record_error(entries: Sequence[tuple[int, int]]) -> None:
        error = [0] * len(columns)
        for coordinate, value in entries:
            error[coordinate] = value
        syndrome = word_syndrome(error, check_rows, q)
        errors_by_syndrome.setdefault(syndrome, set()).add(tuple(error))

    record_error(())
    for weight in range(1, 5):
        for support in combinations(range(len(columns)), weight):
            for values in product(range(1, q), repeat=weight):
                record_error(tuple(zip(support, values)))

    list_size_histogram: dict[int, int] = {}
    maximum_list_size = 0
    coset_count = 0
    for syndrome in product(range(q), repeat=5):
        received = [0] * len(columns)
        basis_values = mat_vec(basis_inverse, syndrome, q)
        for coordinate, value in zip(basis_indices, basis_values):
            received[coordinate] = value
        assert word_syndrome(received, check_rows, q) == syndrome
        decoded_list = tuple(list_decoder(received))
        recovered_errors = {
            tuple(
                (symbol - decoded_symbol) % q
                for symbol, decoded_symbol in zip(received, candidate)
            )
            for candidate in decoded_list
        }
        assert recovered_errors == errors_by_syndrome.get(syndrome, set())
        list_size = len(decoded_list)
        list_size_histogram[list_size] = list_size_histogram.get(list_size, 0) + 1
        maximum_list_size = max(maximum_list_size, list_size)
        coset_count += 1
    error_count = sum(len(errors) for errors in errors_by_syndrome.values())
    return coset_count, maximum_list_size, list_size_histogram, error_count


def exhaustive_unique_decoder(
    columns: Sequence[Vector], received: Sequence[int], q: int, max_errors: int
) -> Vector:
    """Small-parameter oracle used only to regression-test decoder transport."""
    dimension = len(columns[0])
    candidates = []
    for message in product(range(q), repeat=dimension):
        codeword = encode_from_columns(message, columns, q)
        if sum(x != y for x, y in zip(codeword, received)) <= max_errors:
            candidates.append(codeword)
            if len(candidates) > 1:
                break
    if len(candidates) != 1:
        raise ValueError("the received word has no unique codeword in the requested radius")
    return candidates[0]


def canonical_projective_key(points: Iterable[Sequence[int]], q: int) -> tuple[Vector, ...]:
    """Return a canonical key for an unordered projective point set.

    The minimum is taken over all ordered projective frames contained in the
    set.  Consequently, two simple spanning point sets have the same key if
    and only if they are projectively equivalent.
    """
    pts = tuple(normalize_projective(v, q) for v in points)
    if len(set(pts)) != len(pts):
        raise ValueError("the projective point set contains repetitions")
    r = len(pts[0])
    best: tuple[Vector, ...] | None = None
    for indices in permutations(range(len(pts)), r + 1):
        frame = tuple(pts[i] for i in indices)
        if not is_frame(frame, q):
            continue
        T = frame_to_standard_projectivity(frame, q)
        image = tuple(sorted(normalize_projective(mat_vec(T, p, q), q) for p in pts))
        if best is None or image < best:
            best = image
    if best is None:
        raise ValueError("the point set contains no projective frame")
    return best


def standard_cremona(point: Sequence[int], q: int) -> Vector:
    """Standard Cremona involution on the projective torus."""
    p = tuple(x % q for x in point)
    if not all(p):
        raise ValueError("the torus chart requires all coordinates to be nonzero")
    total = 1
    for x in p:
        total = total * x % q
    return normalize_projective(tuple(total * inv_mod(x, q) % q for x in p), q)


def line_coordinates(point: Sequence[int], direction: Sequence[int], q: int) -> tuple[int, int] | None:
    """Solve point = s*(1,...,1)+t*direction, if projectively possible."""
    p = tuple(x % q for x in point)
    w = tuple(x % q for x in direction)
    pair = next(((i, j) for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] != w[j]), None)
    if pair is None:
        raise ValueError("the direction point equals the all-one point")
    i, j = pair
    t = (p[i] - p[j]) * inv_mod(w[i] - w[j], q) % q
    s = (p[i] - t * w[i]) % q
    if all((s + t * w[h] - p[h]) % q == 0 for h in range(len(w))):
        return s, t
    return None


def inverse_cremona_line_tangent(
    left: Sequence[int], right: Sequence[int], s: int, t: int, q: int
) -> tuple[Vector, Vector]:
    """Point and tangent direction on the inverse-Cremona image of a line."""
    u = tuple(x % q for x in left)
    w = tuple(x % q for x in right)
    if len(u) != len(w):
        raise ValueError("line endpoints must have the same dimension")
    ell = tuple((s * x + t * y) % q for x, y in zip(u, w))
    r = len(u)
    point = []
    derivative_s = []
    derivative_t = []
    for i in range(r):
        others = [j for j in range(r) if j != i]
        value = 1
        for j in others:
            value = value * ell[j] % q
        point.append(value)
        ds = 0
        dt = 0
        for h in others:
            product = 1
            for j in others:
                if j != h:
                    product = product * ell[j] % q
            ds = (ds + u[h] * product) % q
            dt = (dt + w[h] * product) % q
        derivative_s.append(ds)
        derivative_t.append(dt)
    tangent = tuple(derivative_s if t % q else derivative_t)
    if matrix_rank(columns_to_matrix((tuple(point), tangent)), q) < 2:
        for candidate in (tuple(derivative_s), tuple(derivative_t)):
            if matrix_rank(columns_to_matrix((tuple(point), candidate)), q) == 2:
                tangent = candidate
                break
        else:
            raise ValueError("degenerate inverse-Cremona parametrization")
    return tuple(point), tangent


def inverse_cremona_curve_tangent(
    direction: Sequence[int], s: int, t: int, q: int
) -> tuple[Vector, Vector]:
    """Specialization to a line through the all-one point."""
    return inverse_cremona_line_tangent(tuple([1] * len(direction)), direction, s, t, q)


def line_is_admissible(left: Sequence[int], right: Sequence[int], q: int) -> bool:
    """Whether a line avoids all codimension-two coordinate subspaces."""
    return all(
        (left[i] * right[j] - left[j] * right[i]) % q != 0
        for i in range(len(left))
        for j in range(i + 1, len(left))
    )


def product_linear_forms(forms: Sequence[tuple[int, int]], q: int) -> Vector:
    """Coefficients of a product of forms a*s+b*t, ordered by t-degree."""
    coefficients = [1]
    for a, b in forms:
        updated = [0] * (len(coefficients) + 1)
        for j, value in enumerate(coefficients):
            updated[j] = (updated[j] + a * value) % q
            updated[j + 1] = (updated[j + 1] + b * value) % q
        coefficients = updated
    return tuple(coefficients)


def inverse_cremona_curve_matrix(left: Sequence[int], right: Sequence[int], q: int) -> Matrix:
    """Coefficient matrix D with d(s,t)=D*(s^k,...,t^k)^T."""
    r = len(left)
    return [
        list(product_linear_forms(tuple((left[j], right[j]) for j in range(r) if j != i), q))
        for i in range(r)
    ]


def symmetric_power_matrix(A: Matrix, k: int, q: int) -> Matrix:
    """Ambient action induced by a 2x2 projectivity on the degree-k Veronese curve."""
    if len(A) != 2 or any(len(row) != 2 for row in A):
        raise ValueError("parameter projectivity must be 2x2")
    first, second = tuple(A[0]), tuple(A[1])
    rows = []
    for j in range(k + 1):
        forms = tuple([first] * (k - j) + [second] * j)
        rows.append(list(product_linear_forms(forms, q)))
    return rows


def parameter_map_to_infinity(s: int, t: int, q: int) -> Matrix:
    """Return A in PGL(2,q) carrying [s:t] to [0:1]."""
    s, t = s % q, t % q
    if s == 0:
        if t == 0:
            raise ValueError("zero parameter vector")
        return [[1, 0], [0, 1]]
    return [[t, -s % q], [inv_mod(s, q), 0]]


def boundary_global_rl_tangent_pairs(
    S: Sequence[int], k: int, lam: int, q: int
) -> list[dict[str, object]]:
    """Global RL tangent_pairs at n=k+3 for a normalized g_(k+1)=0 MDS point set.

    The target evaluation set is unrestricted.  An admissible tangent pair
    consists of an
    excluded source evaluation point, an alternative rational normal curve
    through the other k+4 columns, and a point of tangency containing the
    excluded column.
    """
    if len(S) != k + 3:
        raise ValueError("the boundary test requires n=k+3")
    source = theorem22_case_i_points(S, k, lam, q)
    if not is_arc(source, q):
        raise ValueError("the torus-chart test currently requires an MDS source")
    p_zero, p_lam = source[-2:]
    evaluations = source[:-2]
    tangent_pairs: list[dict[str, object]] = []
    for excluded in range(len(S)):
        retained_values = tuple(a for i, a in enumerate(S) if i != excluded)
        frame = tuple(p for i, p in enumerate(evaluations) if i != excluded)
        H = frame_to_standard_projectivity(frame, q)
        z_zero = normalize_projective(mat_vec(H, p_zero, q), q)
        z_lam = normalize_projective(mat_vec(H, p_lam, q), q)
        z_excluded = normalize_projective(mat_vec(H, evaluations[excluded], q), q)
        basis_values, t_value = retained_values[: k + 1], retained_values[k + 1]
        r_value = sum(basis_values) % q
        predicted_zero = normalize_projective(
            tuple((b - r_value) * (t_value - b) % q for b in basis_values), q
        )
        predicted_lam = normalize_projective(
            tuple((b - r_value + lam) * (t_value - b) % q for b in basis_values), q
        )
        assert z_zero == predicted_zero
        assert z_lam == predicted_lam
        c_zero = standard_cremona(z_zero, q)
        c_lam = standard_cremona(z_lam, q)
        coordinates_lam = line_coordinates(c_lam, c_zero, q)
        if coordinates_lam is None:
            continue
        if len(set(c_zero)) != k + 1:
            continue
        tangent_parameters: list[tuple[str, int, int]] = [("frame_sum", 1, 0), ("P_0", 0, 1)]
        tangent_parameters.append(("P_lambda", coordinates_lam[0], coordinates_lam[1]))
        tangent_parameters.extend((f"frame_e_{i}", -c_zero[i] % q, 1) for i in range(k + 1))
        for label, s, t in tangent_parameters:
            curve_point, tangent = inverse_cremona_curve_tangent(c_zero, s, t, q)
            if matrix_rank(columns_to_matrix((curve_point, tangent, z_excluded)), q) <= 2:
                tangent_pairs.append(
                    {
                        "excluded_index": excluded,
                        "excluded_value": S[excluded],
                        "tangency_point": label,
                        "line_direction": c_zero,
                    }
                )
    return tangent_pairs


def minimal_barycentric_tangent_condition(
    u: Sequence[int],
    v: Sequence[int],
    z: Sequence[int],
    tangent_position: int,
    q: int,
) -> bool:
    """Closed arithmetic form of a minimal-length tangent test.

    The coordinates are relative to a basis of k+1 retained points.  The two
    remaining retained points are u and v, while z is the excluded point.
    Positions 0,...,k denote basis points, k+1 denotes u, and k+2 denotes v.
    The arc hypothesis makes every displayed denominator nonzero.
    """
    r = len(u)
    if len(v) != r or len(z) != r:
        raise ValueError("barycentric coordinate vectors must have equal length")
    if not all(x % q for vector in (u, v, z) for x in vector):
        raise ValueError("the barycentric tangent test requires torus coordinates")
    if tangent_position == r:
        return matrix_rank(
            [
                [1] * r,
                [u[i] * inv_mod(v[i], q) % q for i in range(r)],
                [z[i] * inv_mod(u[i], q) % q for i in range(r)],
            ],
            q,
        ) <= 2
    if tangent_position == r + 1:
        return matrix_rank(
            [
                [1] * r,
                [v[i] * inv_mod(u[i], q) % q for i in range(r)],
                [z[i] * inv_mod(v[i], q) % q for i in range(r)],
            ],
            q,
        ) <= 2
    if not 0 <= tangent_position < r:
        raise ValueError("invalid tangent position")
    m = tangent_position
    values = []
    for i in range(r):
        if i == m:
            continue
        numerator = z[i] * (v[m] * u[i] - u[m] * v[i])
        denominator = u[i] * v[i]
        values.append(numerator * inv_mod(denominator, q) % q)
    return len(set(values)) == 1


def k3_evaluation_tangent_condition(
    S: Sequence[int], excluded_value: int, tangent_value: int, lam: int, q: int
) -> bool:
    """Closed two-equation test for an evaluation/evaluation tangent at k=3."""
    if len(S) != 5 or excluded_value == tangent_value:
        raise ValueError("the k=3 test requires two distinct values in a five-set")
    remaining = tuple(a % q for a in S if a not in (excluded_value, tangent_value))
    if len(remaining) != 3:
        raise ValueError("S must contain distinct evaluation values")
    m, t = tangent_value % q, excluded_value % q
    s1 = sum(remaining) % q
    s2 = sum(remaining[i] * remaining[j] for i in range(3) for j in range(i + 1, 3)) % q
    s3 = remaining[0] * remaining[1] * remaining[2] % q
    linear_condition = (lam - m - sum(S)) % q
    cubic_condition = (
        m**3 + m * m * s1 - m * s1 * t + m * s2 - s1 * t * t - s3
    ) % q
    return linear_condition == 0 and cubic_condition == 0


def k3_minimal_closed_excluded_labels(S: Sequence[int], lam: int, q: int) -> tuple[str, ...]:
    """Complete closed classification of minimal global RL tangent_pairs at k=3.

    If an admissible tangent pair exists, its tangency value is forced to be
    m=lambda-sum(S).  The returned labels are precisely the columns that may
    be excluded while the other six columns lie on the corresponding twisted
    cubic and tangent.
    """
    if len(S) != 5 or len(set(S)) != 5:
        raise ValueError("the k=3 closed classification requires a five-set")
    total_1 = sum(S) % q
    tangent_value = (lam - total_1) % q
    if tangent_value not in S:
        return ()
    total_2 = sum(S[i] * S[j] for i in range(5) for j in range(i + 1, 5)) % q
    total_3 = sum(
        S[i] * S[j] * S[h]
        for i in range(5)
        for j in range(i + 1, 5)
        for h in range(j + 1, 5)
    ) % q
    m = tangent_value
    labels: list[str] = []
    phi_p0 = (total_3 + 3 * total_1 * m * m + 2 * m**3) % q
    if phi_p0 == 0:
        labels.append("P_0")
    phi_plambda = (
        total_3 - total_2 * (total_1 + m) + total_1 * total_1 * m - 2 * m**3
    ) % q
    if phi_plambda == 0:
        labels.append("P_lambda")
    for t in S:
        if t == m:
            continue
        phi_evaluation = (
            2 * t**3
            + (4 * m - 2 * total_1) * t * t
            + (total_2 - 3 * total_1 * m + 2 * m * m) * t
            + 2 * total_2 * m
            - total_3
            - total_1 * m * m
            + 2 * m**3
        ) % q
        if phi_evaluation == 0:
            labels.append(f"evaluation_{t}")
    return tuple(labels)


def k3_unified_candidate_polynomial(
    S: Sequence[int], lam: int, q: int
) -> tuple[tuple[tuple[str, int], ...], tuple[int, int, int, int]]:
    """Return the six unified tangent-pair candidates and their cubic.

    For m=lambda-e_1(S) in S, the two special-column residuals are the
    negatives of the evaluation residual at -2m and e_1(S)-m.  Their zero
    conditions are unchanged.  If C is the resulting six-set, the returned
    coefficients define the symmetric cubic

        2 X^3 - c_1 X^2 + (c_2-c_1^2/4) X
        + c_1(c_2-c_1^2/4)/2 - c_3.

    Labels record which excluded column corresponds to each parameter.  An
    empty candidate tuple means that the forced tangent value is not in S.
    """
    if len(S) != 5 or len(set(S)) != 5:
        raise ValueError("the k=3 unified criterion requires a five-set")
    if q % 2 == 0:
        raise ValueError("the symmetric cubic requires odd characteristic")
    total_1 = sum(S) % q
    m = (lam - total_1) % q
    if m not in S:
        return (), (0, 0, 0, 0)
    candidates = tuple(
        [("P_0", (-2 * m) % q), ("P_lambda", (total_1 - m) % q)]
        + [(f"evaluation_{t}", t % q) for t in S if t != m]
    )
    coefficients = symmetric_six_set_cubic(tuple(value for _, value in candidates), q)
    return candidates, coefficients


def symmetric_six_set_cubic(values: Sequence[int], q: int) -> tuple[int, int, int, int]:
    """Return the coefficient tuple of the symmetric cubic Psi_C."""
    if len(values) != 6 or len(set(values)) != 6:
        raise ValueError("the symmetric cubic requires a six-set")
    if q % 2 == 0:
        raise ValueError("the symmetric cubic requires odd characteristic")
    c1 = sum(values) % q
    c2 = sum(values[i] * values[j] for i in range(6) for j in range(i + 1, 6)) % q
    c3 = sum(
        values[i] * values[j] * values[h]
        for i in range(6)
        for j in range(i + 1, 6)
        for h in range(j + 1, 6)
    ) % q
    quarter = inv_mod(4, q)
    half = inv_mod(2, q)
    linear = (c2 - c1 * c1 * quarter) % q
    return (
        2 % q,
        (-c1) % q,
        linear,
        (c1 * linear * half - c3) % q,
    )


def evaluate_cubic(coefficients: Sequence[int], value: int, q: int) -> int:
    """Evaluate a cubic coefficient tuple by Horner's rule."""
    if len(coefficients) != 4:
        raise ValueError("a cubic needs four coefficients")
    result = 0
    for coefficient in coefficients:
        result = (result * value + coefficient) % q
    return result


def has_balanced_triple_partition(values: Sequence[int], q: int) -> bool:
    """Test whether a six-set splits into two triples having the same sum."""
    if len(values) != 6 or len(set(values)) != 6:
        raise ValueError("the balanced-partition test requires a six-set")
    total = sum(values) % q
    return any(
        (2 * sum(values[i] for i in indices) - total) % q == 0
        for indices in combinations(range(6), 3)
    )


def quadratic_character(value: int, q: int) -> int:
    """Quadratic character over an odd prime field."""
    value %= q
    if value == 0:
        return 0
    symbol = pow(value, (q - 1) // 2, q)
    if symbol == 1:
        return 1
    if symbol == q - 1:
        return -1
    raise ValueError("q must be an odd prime")


def k3_closed_rooted_translation_count(q: int) -> int:
    """Closed number of rooted six-set translation classes in characteristic >5."""
    if q <= 5 or any(q % divisor == 0 for divisor in range(2, int(q**0.5) + 1)):
        raise ValueError("the closed orbit count requires a prime characteristic greater than 5")
    open_projective_count = (
        q**3
        - 19 * q**2
        + 121 * q
        - 269
        - 30 * quadratic_character(-1, q)
        - 20 * quadratic_character(-3, q)
    )
    numerator = (q - 1) * open_projective_count
    assert numerator % 120 == 0
    return numerator // 120


def k3_normalized_double_root_count(q: int) -> int:
    """Count normalized six-sets whose candidate cubic has roots 0 and 1.

    The six-set is ``{0, 1} union A``.  The returned number counts four-subsets
    ``A`` of the remaining field elements, not ordered four-tuples.
    """
    if q <= 5 or any(q % divisor == 0 for divisor in range(2, int(q**0.5) + 1)):
        raise ValueError("the normalized double-root count requires a prime greater than 5")
    count = 0
    for values in combinations(range(2, q), 4):
        a1 = sum(values) % q
        a2 = sum(values[i] * values[j] for i in range(4) for j in range(i + 1, 4)) % q
        a3 = sum(
            values[i] * values[j] * values[h]
            for i in range(4)
            for j in range(i + 1, 4)
            for h in range(j + 1, 4)
        ) % q
        if (4 * a2 - (a1 - 1) * (a1 + 3)) % q != 0:
            continue
        if (4 * a3 - (a1 - 1) ** 2) % q == 0:
            count += 1
    return count


def k3_cm_minus8_trace(q: int) -> int:
    """Frobenius trace of ``y^2 = x^3 - 4*x^2 + 2*x`` over a prime field."""
    if q <= 5 or any(q % divisor == 0 for divisor in range(2, int(q**0.5) + 1)):
        raise ValueError("the CM trace requires a prime characteristic greater than 5")
    return -sum(
        quadratic_character(x**3 - 4 * x**2 + 2 * x, q)
        for x in range(q)
    )


def k3_cm_minus8_trace_square(q: int) -> int:
    """Closed CM formula for the square of the Frobenius trace.

    For split primes, ``q = c^2 + 2*d^2`` and the answer is ``4*c^2``.
    For inert primes the trace is zero.  This avoids enumerating the points of
    the elliptic curve when only the trace square is needed.
    """
    if q <= 5 or any(q % divisor == 0 for divisor in range(2, isqrt(q) + 1)):
        raise ValueError("the CM trace formula requires a prime greater than 5")
    return k3_cm_minus8_trace_square_prime_power(q, 1)


def k3_cm_minus8_trace_square_prime_power(p: int, r: int) -> int:
    """CM trace square over ``F_(p^r)`` using only integer arithmetic."""
    if p <= 5 or any(p % divisor == 0 for divisor in range(2, isqrt(p) + 1)):
        raise ValueError("p must be a prime greater than 5")
    if r < 1:
        raise ValueError("r must be positive")
    if p % 8 in (5, 7):
        return 0 if r % 2 else 4 * p**r
    c = None
    for d in range(1, isqrt(p // 2) + 1):
        c_squared = p - 2 * d * d
        if c_squared <= 0:
            break
        candidate = isqrt(c_squared)
        if candidate * candidate == c_squared:
            c = candidate
            break
    if c is None:
        raise AssertionError("a split prime must be represented by c^2 + 2*d^2")
    previous, current = 2, 2 * c
    for _ in range(2, r + 1):
        previous, current = current, 2 * c * current - p * previous
    return current**2


def k3_closed_normalized_double_root_count_prime_power(p: int, r: int) -> int:
    """Closed normalized overlap count over ``F_(p^r)``."""
    trace_square = k3_cm_minus8_trace_square_prime_power(p, r)
    q = p**r
    chi = lambda value: quadratic_character(value, p) ** r
    numerator = (
        q**2
        - (16 + chi(-2)) * q
        + 61
        + 6 * chi(-1)
        + 9 * chi(2)
        + 8 * chi(-2)
        + trace_square
    )
    assert numerator % 24 == 0
    return numerator // 24


def k3_closed_normalized_double_root_count(q: int) -> int:
    """Closed formula for the normalized double-tangent-pair overlap ``M_q``.

    The non-polynomial term is the square of the Frobenius trace of the
    discriminant ``-8`` CM elliptic curve used by the K3 character sum.
    """
    return k3_closed_normalized_double_root_count_prime_power(q, 1)


def k3_closed_global_rl_pair_count(q: int) -> int:
    """Closed global RL-pair count, including the CM K3 overlap term."""
    tangent_pair_count = 30 * k3_closed_rooted_translation_count(q)
    double_tangent_presentations = (
        15 * (q - 1) * k3_closed_normalized_double_root_count(q)
    )
    return tangent_pair_count - double_tangent_presentations


def k3_closed_involution_fixed_rl_pair_count(q: int) -> int:
    """Closed RL-pair count fixed by the special-column involution."""
    numerator = q - 6 - quadratic_character(-1, q) - 2 * quadratic_character(2, q)
    assert numerator % 8 == 0
    return 3 * (q - 1) * (numerator // 8)


def k3_closed_monomial_subgroup_orbit_count(q: int) -> int:
    """Burnside count for the scaling and special-involution subgroup."""
    numerator = (
        k3_closed_global_rl_pair_count(q)
        + k3_closed_involution_fixed_rl_pair_count(q)
    )
    denominator = 2 * (q - 1)
    assert numerator % denominator == 0
    return numerator // denominator


def k3_case_i_is_mds_by_triple_sums(
    S: Sequence[int], lam: int, q: int
) -> bool:
    """Closed MDS test for the normalized g_(k+1)=0 family at k=3."""
    if len(S) < 3 or len(set(S)) != len(S) or lam % q == 0:
        return False
    triple_sums = {
        sum(triple) % q for triple in combinations(tuple(S), 3)
    }
    return 0 not in triple_sums and lam % q not in triple_sums


def k3_case_i_is_mds(S: Sequence[int], lam: int, q: int) -> bool:
    """Closed MDS test for the seven-column minimal family at k=3."""
    return len(S) == 5 and k3_case_i_is_mds_by_triple_sums(S, lam, q)


def k3_boundary_case_i_is_mds(S: Sequence[int], lam: int, q: int) -> bool:
    """Closed MDS test for the eight-column boundary family at k=3."""
    return len(S) == 6 and k3_case_i_is_mds_by_triple_sums(S, lam, q)


def boundary_additive_signature(
    S: Sequence[int], lam: int, q: int
) -> tuple[int, int, int]:
    """Return the marked additive signature ``(eta,tau,mu)``.

    Here ``eta`` records whether ``sum(S)=lambda``, ``tau`` is the size of
    the restricted three-fold sumset, and ``mu`` is the largest multiplicity
    of an unordered two-sum.  All three entries are invariant under the
    genuine scaling action ``(S,lambda)->(uS,u lambda)``.
    """
    evaluation_set = tuple(value % q for value in S)
    if len(set(evaluation_set)) != len(evaluation_set) or lam % q == 0:
        raise ValueError("the additive signature requires distinct points and nonzero lambda")
    eta = int(sum(evaluation_set) % q == lam % q)
    restricted_triple_sums = {
        sum(triple) % q for triple in combinations(evaluation_set, 3)
    }
    pair_multiplicities = {value: 0 for value in range(q)}
    for left, right in combinations(evaluation_set, 2):
        pair_multiplicities[(left + right) % q] += 1
    return eta, len(restricted_triple_sums), max(pair_multiplicities.values())


def minimal_global_rl_tangent_pairs(
    S: Sequence[int], k: int, lam: int, q: int
) -> list[dict[str, object]]:
    """Classify global RL equivalence at the minimal length n=k+2.

    For each possible excluded column, the remaining k+3 MDS columns determine
    a unique rational normal curve.  The excluded column gives an RL code
    exactly when it lies on a tangent at one of those k+3 curve points.
    """
    if len(S) != k + 2:
        raise ValueError("the minimal-length test requires n=k+2")
    source = theorem22_case_i_points(S, k, lam, q)
    if not is_arc(source, q):
        raise ValueError("the minimal global criterion currently requires an MDS source")
    labels = tuple([f"evaluation_{a}" for a in S] + ["P_0", "P_lambda"])
    tangent_pairs: list[dict[str, object]] = []
    r = k + 1
    for excluded in range(len(source)):
        included_indices = tuple(i for i in range(len(source)) if i != excluded)
        included = tuple(source[i] for i in included_indices)
        basis = included[:r]
        H = mat_inv(columns_to_matrix(basis), q)
        transformed = tuple(normalize_projective(mat_vec(H, p, q), q) for p in included)
        z_excluded = normalize_projective(mat_vec(H, source[excluded], q), q)
        c_left = standard_cremona(transformed[r], q)
        c_right = standard_cremona(transformed[r + 1], q)
        if not line_is_admissible(c_left, c_right, q):
            continue
        tangent_parameters: list[tuple[int, str, int, int]] = [
            (included_indices[r], labels[included_indices[r]], 1, 0),
            (included_indices[r + 1], labels[included_indices[r + 1]], 0, 1),
        ]
        tangent_parameters.extend(
            (included_indices[i], labels[included_indices[i]], -c_right[i] % q, c_left[i])
            for i in range(r)
        )
        for tangent_position, (tangent_index, tangent_label, s, t) in enumerate(
            tangent_parameters[2:] + tangent_parameters[:2]
        ):
            # Reorder from (u,v,basis points) to (basis points,u,v) so that
            # tangent_position agrees with the closed barycentric criterion.
            curve_point, tangent = inverse_cremona_line_tangent(c_left, c_right, s, t, q)
            geometric_condition = (
                matrix_rank(columns_to_matrix((curve_point, tangent, z_excluded)), q) <= 2
            )
            arithmetic_condition = minimal_barycentric_tangent_condition(
                transformed[r], transformed[r + 1], z_excluded, tangent_position, q
            )
            assert geometric_condition == arithmetic_condition
            if k == 3 and excluded < len(S):
                if tangent_index < len(S):
                    closed_k3_condition = k3_evaluation_tangent_condition(
                        S, S[excluded], S[tangent_index], lam, q
                    )
                else:
                    # The two special-point tangencies are ruled out by the
                    # degree-four interpolation argument in the paper.
                    closed_k3_condition = False
                assert geometric_condition == closed_k3_condition
            if geometric_condition:
                curve_matrix = inverse_cremona_curve_matrix(c_left, c_right, q)
                J = mat_inv(curve_matrix, q)
                parameter_map = parameter_map_to_infinity(s, t, q)
                K = symmetric_power_matrix(parameter_map, k, q)
                projectivity = mat_mul(mat_mul(K, J, q), H, q)
                images = tuple(
                    normalize_projective(mat_vec(projectivity, point, q), q) for point in source
                )
                infinity = tuple([0] * k + [1])
                assert images[tangent_index] == infinity
                target_S = []
                for index in included_indices:
                    if index == tangent_index:
                        continue
                    image = images[index]
                    if image[0] == 0:
                        raise AssertionError("a non-tangency curve point mapped to infinity")
                    value = image[1] * inv_mod(image[0], q) % q
                    assert image == normalize_projective(nu(value, k, q), q)
                    target_S.append(value)
                off_point = images[excluded]
                assert all(x == 0 for x in off_point[: k - 1])
                assert off_point[k - 1] != 0
                delta = off_point[k] * inv_mod(off_point[k - 1], q) % q
                target_S_tuple = tuple(sorted(target_S))
                assert len(target_S_tuple) == k + 2
                assert len(set(target_S_tuple)) == k + 2
                target = roth_lempel_points(target_S_tuple, k, delta, q)
                assert set(images) == {normalize_projective(point, q) for point in target}
                monomial_map = monomial_map_from_projectivity(
                    source, target, projectivity, q
                )
                tangent_pairs.append(
                    {
                        "excluded_index": excluded,
                        "excluded_label": labels[excluded],
                        "tangency_point": tangent_label,
                        "tangent_parameter": (s, t),
                        "line_left": c_left,
                        "line_right": c_right,
                        "target_S": target_S_tuple,
                        "delta": delta,
                        "projectivity": projectivity,
                        "monomial_map": monomial_map,
                    }
                )
    if k == 3:
        tangent_label = f"evaluation_{(lam - sum(S)) % q}"
        closed_pairs = {
            (excluded_label, tangent_label)
            for excluded_label in k3_minimal_closed_excluded_labels(S, lam, q)
        }
        geometric_pairs = {
            (str(pair["excluded_label"]), str(pair["tangency_point"]))
            for pair in tangent_pairs
        }
        assert geometric_pairs == closed_pairs
    return tangent_pairs


def nu(a: int, k: int, q: int) -> Vector:
    return tuple(pow(a, j, q) for j in range(k + 1))


def theorem22_case_i_points(S: Sequence[int], k: int, lam: int, q: int) -> tuple[Vector, ...]:
    """Columns of the normalized matrix in Remark 6(2)."""
    e_k_minus_1 = tuple([0] * (k - 1) + [1, 0])
    second = tuple([0] * (k - 1) + [1, lam % q])
    return tuple(nu(a, k, q) for a in S) + (e_k_minus_1, second)


def case_i_structured_parity_check(
    S: Sequence[int], k: int, lam: int, q: int
) -> tuple[Vector, ...]:
    """Explicit parity check for the normalized g_(k+1)=0 family.

    Put R=|S|-k+1.  With w_a=prod_{b!=a}(a-b)^(-1) and e1=sum(S),
    the evaluation columns are w_a*(1,a,...,a^(R-1)).  The two extension
    columns have R-2 leading zeros followed respectively by
    lambda^(-1)*(1,e1-lambda) and -lambda^(-1)*(1,e1).
    """
    evaluation_set = tuple(value % q for value in S)
    redundancy = len(evaluation_set) - k + 1
    if redundancy < 3 or len(set(evaluation_set)) != len(evaluation_set):
        raise ValueError(
            "the structured formula requires at least k+2 distinct evaluation points"
        )
    if lam % q == 0:
        raise ValueError("lambda must be nonzero")
    inverse_lambda = inv_mod(lam, q)
    first_sum = sum(evaluation_set) % q
    check_columns = []
    for a in evaluation_set:
        derivative_value = 1
        for b in evaluation_set:
            if b != a:
                derivative_value = derivative_value * (a - b) % q
        barycentric_weight = inv_mod(derivative_value, q)
        check_columns.append(
            tuple(
                barycentric_weight * pow(a, exponent, q) % q
                for exponent in range(redundancy)
            )
        )
    check_columns.append(
        (0,) * (redundancy - 2)
        + (inverse_lambda, (first_sum - lam) * inverse_lambda % q)
    )
    check_columns.append(
        (0,) * (redundancy - 2)
        + (-inverse_lambda % q, -first_sum * inverse_lambda % q)
    )
    parity_check = tuple(tuple(row) for row in zip(*check_columns))
    generator = columns_to_matrix(theorem22_case_i_points(S, k, lam, q))
    assert all(
        sum(generator[row][column] * parity_check[check][column]
            for column in range(len(check_columns))) % q
        == 0
        for row in range(k + 1)
        for check in range(redundancy)
    )
    assert matrix_rank([list(row) for row in parity_check], q) == redundancy
    return parity_check


def boundary_case_i_structured_parity_check(
    S: Sequence[int], k: int, lam: int, q: int
) -> tuple[Vector, ...]:
    """Four-row specialization of :func:`case_i_structured_parity_check`."""
    if len(S) != k + 3:
        raise ValueError("the boundary formula requires k+3 evaluation points")
    return case_i_structured_parity_check(S, k, lam, q)


def canonical_case_i_scaling_pair(
    S: Sequence[int], lam: int, q: int
) -> tuple[tuple[int, ...], int]:
    """Canonicalize ``(S,lambda)`` under ``(S,lambda)->(uS,u lambda)``.

    This is a genuine monomial equivalence of the normalized g_(k+1)=0
    family.  In contrast, translations need not preserve its marked tangent
    pair and are deliberately excluded.
    """
    return min(
        (tuple(sorted(u * a % q for a in S)), u * lam % q)
        for u in range(1, q)
    )


def radius_three_syndrome_histogram(
    parity_check: Sequence[Sequence[int]], q: int
) -> tuple[int, int, int, dict[int, int]]:
    """Count exact weight-at-most-three lists in every redundancy-four coset."""
    check_rows = tuple(tuple(entry % q for entry in row) for row in parity_check)
    if len(check_rows) != 4 or not check_rows:
        raise ValueError("the histogram requires a four-row parity check")
    length = len(check_rows[0])
    if any(len(row) != length for row in check_rows):
        raise ValueError("parity-check rows have inconsistent lengths")
    check_columns = tuple(
        tuple(check_rows[row][coordinate] for row in range(4))
        for coordinate in range(length)
    )
    syndrome_counts: dict[Vector, int] = {(0, 0, 0, 0): 1}

    def record(support: Sequence[int], values: Sequence[int]) -> None:
        syndrome = tuple(
            sum(value * check_columns[coordinate][row]
                for coordinate, value in zip(support, values)) % q
            for row in range(4)
        )
        syndrome_counts[syndrome] = syndrome_counts.get(syndrome, 0) + 1

    for weight in range(1, 4):
        for support in combinations(range(length), weight):
            for values in product(range(1, q), repeat=weight):
                record(support, values)

    total_cosets = q**4
    histogram: dict[int, int] = {0: total_cosets - len(syndrome_counts)}
    for count in syndrome_counts.values():
        histogram[count] = histogram.get(count, 0) + 1
    if histogram[0] == 0:
        del histogram[0]
    error_count = sum(size * frequency for size, frequency in histogram.items())
    maximum_list_size = max(histogram)
    return total_cosets, error_count, maximum_list_size, dict(sorted(histogram.items()))


def projective_points_prime_field(dimension: int, q: int) -> tuple[Vector, ...]:
    """Enumerate normalized points of ``PG(dimension-1,q)`` over a prime field."""
    if dimension < 1:
        raise ValueError("the vector-space dimension must be positive")
    points = []
    for pivot in range(dimension):
        prefix = (0,) * pivot + (1,)
        for suffix in product(range(q), repeat=dimension - pivot - 1):
            points.append(prefix + suffix)
    expected = sum(q**power for power in range(dimension))
    assert len(points) == expected
    return tuple(points)


def codimension_four_projective_hole_count(
    parity_check: Sequence[Sequence[int]],
    q: int,
    *,
    projective_points: Sequence[Vector] | None = None,
) -> int:
    """Count projective syndromes outside every three-column check plane.

    NumPy is imported only for this exhaustive finite-field audit.  The
    decoder and all core equivalence routines remain dependency-free.
    """
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("the projective-hole audit requires NumPy") from exc

    check_rows = tuple(tuple(entry % q for entry in row) for row in parity_check)
    if len(check_rows) != 4 or any(len(row) != len(check_rows[0]) for row in check_rows):
        raise ValueError("the projective-hole audit requires a four-row parity check")
    length = len(check_rows[0])
    check_columns = tuple(
        tuple(check_rows[row][coordinate] for row in range(4))
        for coordinate in range(length)
    )
    plane_keys = []
    for support in combinations(range(length), 3):
        nullspace = right_nullspace_basis(
            [list(check_columns[coordinate]) for coordinate in support], q
        )
        if len(nullspace) != 1:
            raise ValueError("three check columns do not span a projective plane")
        plane_keys.append(normalize_projective(nullspace[0], q))
    if len(set(plane_keys)) != len(plane_keys):
        raise ValueError("distinct coordinate triples produced the same check plane")

    points = (
        projective_points_prime_field(4, q)
        if projective_points is None
        else tuple(tuple(entry % q for entry in point) for point in projective_points)
    )
    point_matrix = np.asarray(points, dtype=np.int64)
    plane_matrix = np.asarray(plane_keys, dtype=np.int64)
    covered = np.any((point_matrix @ plane_matrix.T) % q == 0, axis=1)
    return int(np.count_nonzero(~covered))


def codimension_five_projective_hole_count(
    parity_check: Sequence[Sequence[int]],
    q: int,
    *,
    projective_points: Sequence[Vector] | None = None,
) -> int:
    """Count projective syndromes outside every four-column check hyperplane."""
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("the projective-hole audit requires NumPy") from exc

    check_rows = tuple(tuple(entry % q for entry in row) for row in parity_check)
    if len(check_rows) != 5 or any(len(row) != len(check_rows[0]) for row in check_rows):
        raise ValueError("the projective-hole audit requires a five-row parity check")
    length = len(check_rows[0])
    check_columns = tuple(
        tuple(check_rows[row][coordinate] for row in range(5))
        for coordinate in range(length)
    )
    hyperplane_keys = codimension_five_check_hyperplane_keys(check_columns, q)
    if len(set(hyperplane_keys)) != len(hyperplane_keys):
        raise ValueError("distinct coordinate quadruples produced the same check hyperplane")

    points = (
        projective_points_prime_field(5, q)
        if projective_points is None
        else tuple(tuple(entry % q for entry in point) for point in projective_points)
    )
    point_matrix = np.asarray(points, dtype=np.int64)
    hyperplane_matrix = np.asarray(hyperplane_keys, dtype=np.int64)
    covered = np.any((point_matrix @ hyperplane_matrix.T) % q == 0, axis=1)
    return int(np.count_nonzero(~covered))


def codimension_five_check_hyperplane_keys(
    check_columns: Sequence[Vector], q: int
) -> tuple[Vector, ...]:
    """Return the hyperplanes spanned by all check-column quadruples."""
    if not check_columns or any(len(column) != 5 for column in check_columns):
        raise ValueError("five-dimensional check columns are required")
    hyperplane_keys = []
    for support in combinations(range(len(check_columns)), 4):
        nullspace = right_nullspace_basis(
            [list(check_columns[coordinate]) for coordinate in support], q
        )
        if len(nullspace) != 1:
            raise ValueError("four check columns do not span a projective hyperplane")
        hyperplane_keys.append(normalize_projective(nullspace[0], q))
    return tuple(hyperplane_keys)


def polynomial_coefficients_from_roots(
    roots: Sequence[int], q: int
) -> Vector:
    """Return ascending coefficients of the monic polynomial with given roots."""
    coefficients = [1]
    for root in roots:
        updated = [0] * (len(coefficients) + 1)
        for degree, coefficient in enumerate(coefficients):
            updated[degree] = (updated[degree] - root * coefficient) % q
            updated[degree + 1] = (updated[degree + 1] + coefficient) % q
        coefficients = updated
    return tuple(coefficients)


def distance_six_completion_hyperplane_keys(
    S: Sequence[int], lam: int, q: int
) -> tuple[Vector, ...]:
    """Factorized completion arrangement for the k=3, redundancy-five family."""
    evaluation_set = tuple(value % q for value in S)
    if len(evaluation_set) != 7 or len(set(evaluation_set)) != 7:
        raise ValueError("seven distinct evaluation points are required")
    if lam % q == 0:
        raise ValueError("lambda must be nonzero")
    first_sum = sum(evaluation_set) % q
    polynomial_normals: list[Vector] = []
    for support in combinations(evaluation_set, 4):
        polynomial_normals.append(polynomial_coefficients_from_roots(support, q))
    for support in combinations(evaluation_set, 3):
        support_sum = sum(support) % q
        for tangent_parameter in (first_sum - lam, first_sum):
            fourth_root = (tangent_parameter - support_sum) % q
            polynomial_normals.append(
                polynomial_coefficients_from_roots(support + (fourth_root,), q)
            )
    for support in combinations(evaluation_set, 2):
        quadratic = polynomial_coefficients_from_roots(support, q)
        polynomial_normals.append(quadratic + (0, 0))
    normalized = tuple(normalize_projective(normal, q) for normal in polynomial_normals)
    assert len(normalized) == 126
    return normalized


def completion_arrangement_first_coefficients(
    normals: Sequence[Vector], q: int
) -> tuple[int, int, dict[int, int], int]:
    """Return b2, b3, the rank-two multiplicity histogram, and rank-three flats.

    The characteristic polynomial begins
    ``T^5-len(normals)*T^4+b2*T^3-b3*T^2``.
    """
    normal_count = len(normals)
    if not normals or any(len(normal) != 5 for normal in normals):
        raise ValueError("five-dimensional arrangement normals are required")

    pair_keys: dict[tuple[int, int], tuple[Vector, ...]] = {}
    pair_counts: dict[tuple[Vector, ...], int] = {}
    for left, right in combinations(range(normal_count), 2):
        key = row_space_key((normals[left], normals[right]), q)
        if len(key) != 2:
            raise ValueError("arrangement hyperplanes must be distinct")
        pair_keys[left, right] = key
        pair_counts[key] = pair_counts.get(key, 0) + 1

    pair_mobius: dict[tuple[Vector, ...], int] = {}
    rank_two_multiplicities: dict[int, int] = {}
    for key, pair_count in pair_counts.items():
        multiplicity = (1 + isqrt(1 + 8 * pair_count)) // 2
        if multiplicity * (multiplicity - 1) // 2 != pair_count:
            raise AssertionError("rank-two pair count is not triangular")
        pair_mobius[key] = multiplicity - 1
        rank_two_multiplicities[multiplicity] = (
            rank_two_multiplicities.get(multiplicity, 0) + 1
        )
    b2 = sum(pair_mobius.values())

    rank_three_masks: dict[tuple[Vector, ...], int] = {}
    for first, second, third in combinations(range(normal_count), 3):
        key = row_space_key(
            (normals[first], normals[second], normals[third]), q
        )
        if len(key) != 3:
            continue
        rank_three_masks[key] = (
            rank_three_masks.get(key, 0)
            | (1 << first)
            | (1 << second)
            | (1 << third)
        )

    b3 = 0
    for mask in rank_three_masks.values():
        indices = []
        while mask:
            lowest_bit = mask & -mask
            indices.append(lowest_bit.bit_length() - 1)
            mask -= lowest_bit
        contained_rank_two_flats = {
            pair_keys[min(left, right), max(left, right)]
            for left, right in combinations(indices, 2)
        }
        mobius_value = -(
            1
            - len(indices)
            + sum(pair_mobius[key] for key in contained_rank_two_flats)
        )
        b3 -= mobius_value
    return b2, b3, dict(sorted(rank_two_multiplicities.items())), len(rank_three_masks)


def theorem22_case_ii_points(
    S: Sequence[int], k: int, alpha: int, beta: int, q: int
) -> tuple[Vector, ...]:
    """Normalized columns for g_{k+1} != 0.

    alpha = g_{k-1}/g_{k+1} and beta = (u_{n+1}v_{n+1}-f_k)/g_{k+1}.
    """
    evaluation = []
    for a in S:
        p = [pow(a, j, q) for j in range(k - 1)]
        p.extend([pow(a, k, q), (pow(a, k + 1, q) + alpha * pow(a, k - 1, q)) % q])
        evaluation.append(tuple(p))
    infinity = tuple([0] * (k - 1) + [1, beta % q])
    extension = tuple([0] * k + [1])
    return tuple(evaluation) + (infinity, extension)


def run_nonzero_branch_threshold_example() -> None:
    """Verify the nonempty MDS example at the degree-gap threshold."""
    q, k = 17, 2
    S = (0, 1, 2, 5, 7, 11, 14)
    alpha, beta = 12, 6
    points = theorem22_case_ii_points(S, k, alpha, beta, q)
    assert len(points) == len(S) + 2
    assert is_arc(points, q)
    assert len(S) == 2 * k + 3
    print(
        f"q={q}, k={k}, n={len(S)}: nonzero-leading MDS threshold example, "
        f"alpha={alpha}, beta={beta}"
    )


def roth_lempel_points(S: Sequence[int], k: int, delta: int, q: int) -> tuple[Vector, ...]:
    """Columns of RL_{k+1,delta}(S), embedded in PG(k,q)."""
    infinity = tuple([0] * k + [1])
    tangent = tuple([0] * (k - 1) + [1, delta % q])
    return tuple(nu(a, k, q) for a in S) + (infinity, tangent)


def equivalent_rl_deltas(points: Sequence[Vector], S: Sequence[int], k: int, q: int) -> dict[int, Matrix]:
    witnesses: dict[int, Matrix] = {}
    for delta in range(q):
        target = roth_lempel_points(S, k, delta, q)
        witness = projective_equivalence(points, target, q)
        if witness is not None:
            witnesses[delta] = witness
    return witnesses


def equivalent_rl_deltas_canonical(
    points: Sequence[Vector], S: Sequence[int], k: int, q: int
) -> list[int]:
    """Classify all RL parameters using canonical projective-set keys."""
    source_key = canonical_projective_key(points, q)
    matches = []
    for delta in range(q):
        if source_key == canonical_projective_key(roth_lempel_points(S, k, delta, q), q):
            matches.append(delta)
    return matches


def case_i_equivalence_table(S: Sequence[int], k: int, q: int) -> dict[int, list[int]]:
    """Return lambda -> all fixed-S equivalent RL delta values for g_(k+1)=0."""
    rl_by_key: dict[tuple[Vector, ...], list[int]] = {}
    for delta in range(q):
        key = canonical_projective_key(roth_lempel_points(S, k, delta, q), q)
        rl_by_key.setdefault(key, []).append(delta)
    table = {}
    for lam in range(1, q):
        key = canonical_projective_key(theorem22_case_i_points(S, k, lam, q), q)
        table[lam] = rl_by_key.get(key, [])
    return table


def affine_subset_orbit(S: Sequence[int], q: int) -> set[tuple[int, ...]]:
    """Orbit of a prime-field subset under a -> u*a+b, u != 0."""
    return {
        tuple(sorted((u * a + b) % q for a in S))
        for u in range(1, q)
        for b in range(q)
    }


def affine_subset_representatives(q: int, n: int) -> list[tuple[int, ...]]:
    """Return canonical representatives of n-subsets of the prime field."""
    unseen = set(combinations(range(q), n))
    representatives = []
    while unseen:
        representative = min(unseen)
        orbit = affine_subset_orbit(representative, q)
        representatives.append(representative)
        unseen.difference_update(orbit)
    return representatives


def run_source_paper_checks() -> None:
    q, k = 11, 3

    # Remark 6(2): lambda=3 is equivalent exactly for delta in {3,7}.
    S_remark = (1, 2, 5, 6, 9)
    remark = equivalent_rl_deltas(theorem22_case_i_points(S_remark, k, 3, q), S_remark, k, q)
    assert sorted(remark) == [3, 7]
    print("Remark 6(2):", sorted(remark))

    # Example 2: the displayed matrices have lambda=10,8,9, in that order.
    # This order differs from the preceding order of the listed deep holes.
    S_example = (3, 4, 5, 6, 7)
    expected = {10: [], 8: [], 9: [9, 10]}
    for lam in (10, 8, 9):
        witnesses = equivalent_rl_deltas(theorem22_case_i_points(S_example, k, lam, q), S_example, k, q)
        assert sorted(witnesses) == expected[lam]
        print(f"Example 2, lambda={lam}:", sorted(witnesses))

    # Example 3: alpha=g_2/g_4=4 and beta=(u_6v_6-f_3)/g_4=4 in F_11.
    case_ii = theorem22_case_ii_points(S_example, k, alpha=4, beta=4, q=q)
    witnesses = equivalent_rl_deltas(case_ii, S_example, k, q)
    assert not witnesses
    print("Example 3:", sorted(witnesses))


def run_q11_boundary_scan() -> None:
    """Run a fixed-S boundary sample at q=11, k=3, n=k+3.

    Every nonzero lambda and every RL delta are checked for six affine-orbit
    representatives, but the target evaluation set is constrained to equal
    the source S.  This is a regression sample, not a global RL test: the
    marked tangent pair {Q_0,Q_lambda} is not preserved by translations.
    """
    q, k, n = 11, 3, 6
    representatives = affine_subset_representatives(q, n)
    assert len(representatives) == 6
    print(f"q={q}, k={k}, n={n}: {len(representatives)} affine-orbit representatives")
    for S in representatives:
        table = case_i_equivalence_table(S, k, q)
        positives = {lam: deltas for lam, deltas in table.items() if deltas}
        assert not positives
        print(S, positives)


def run_global_mds_boundary_scan(q: int) -> None:
    """Exhaust the global MDS boundary criterion for k=3 and q in {13,17}."""
    if q not in (13, 17):
        raise ValueError("the verified regression counts are available for q=13 or q=17")
    k, n = 3, 6
    mds_count = 0
    positives: list[tuple[tuple[int, ...], int, list[dict[str, object]]]] = []
    for S in combinations(range(q), n):
        for lam in range(1, q):
            points = theorem22_case_i_points(S, k, lam, q)
            if not is_arc(points, q):
                continue
            mds_count += 1
            tangent_pairs = boundary_global_rl_tangent_pairs(S, k, lam, q)
            if tangent_pairs:
                positives.append((S, lam, tangent_pairs))
    expected_mds_counts = {13: 96, 17: 5904}
    assert mds_count == expected_mds_counts[q]
    assert not positives
    print(f"q={q}, k={k}, n={n}: MDS pairs={mds_count}, global RL tangent_pairs={len(positives)}")


def run_boundary_radius_three_decoder() -> None:
    """Verify structured duals and radius-three decoding through distance six."""
    for q, k, redundancy in ((13, 3, 3), (13, 3, 4), (17, 3, 5), (19, 4, 6)):
        S = tuple(range(k + redundancy - 1))
        for lam in range(1, q):
            structured = case_i_structured_parity_check(S, k, lam, q)
            generic = parity_check_from_columns(theorem22_case_i_points(S, k, lam, q), q)
            assert len(structured) == redundancy
            assert (
                matrix_rank([list(row) for row in structured + generic], q)
                == redundancy
            )

    q, k, lam = 17, 3, 1
    distance_six_set = tuple(range(7))
    distance_six_source = theorem22_case_i_points(distance_six_set, k, lam, q)
    assert k3_case_i_is_mds_by_triple_sums(distance_six_set, lam, q)
    assert is_arc(distance_six_source, q)
    distance_six_check = case_i_structured_parity_check(
        distance_six_set, k, lam, q
    )
    distance_six_decoder = prepare_radius_three_syndrome_list_decoder(
        distance_six_source, q, parity_check=distance_six_check
    )
    zero_codeword = tuple(0 for _ in distance_six_source)
    support_tests = 0
    for weight in range(4):
        for support in combinations(range(len(distance_six_source)), weight):
            received = [0] * len(distance_six_source)
            for position, coordinate in enumerate(support):
                received[coordinate] = position + 1
            decoded = distance_six_decoder(received)
            assert zero_codeword in decoded
            support_tests += 1
    print(
        f"distance-six non-RL q={q}, k={k}, S={distance_six_set}, "
        f"lambda={lam}: support-tests={support_tests}, redundancy=5"
    )

    for q, k in ((13, 3), (17, 4), (19, 5)):
        S = tuple(range(k + 3))
        for lam in range(1, q):
            structured = boundary_case_i_structured_parity_check(S, k, lam, q)
            generic = parity_check_from_columns(theorem22_case_i_points(S, k, lam, q), q)
            assert matrix_rank([list(row) for row in structured + generic], q) == 4

    q, k = 13, 3
    S, lam = (0, 1, 2, 3, 4, 5), 1
    source = theorem22_case_i_points(S, k, lam, q)
    assert is_arc(source, q)
    assert not boundary_global_rl_tangent_pairs(S, k, lam, q)
    structured = boundary_case_i_structured_parity_check(S, k, lam, q)
    decoder = prepare_radius_three_syndrome_list_decoder(
        source, q, parity_check=structured
    )
    cosets, maximum_list_size, histogram, error_count = (
        verify_all_radius_three_syndrome_cosets(
            source, decoder, q, parity_check=structured
        )
    )
    expected_error_count = sum(
        len(tuple(combinations(range(len(source)), weight))) * (q - 1) ** weight
        for weight in range(4)
    )
    assert cosets == q**4
    assert error_count == expected_error_count
    print(
        f"boundary non-RL q={q}, k={k}, S={S}, lambda={lam}: "
        f"syndrome-cosets={cosets}, radius-three errors={error_count}, "
        f"max-list={maximum_list_size}, histogram={dict(sorted(histogram.items()))}"
    )


def run_distance_six_radius_four_decoder() -> None:
    """Verify exact radius-four quotient decoding and the first non-RL instance."""
    q = 5
    repetition_columns = tuple((1,) for _ in range(6))
    repetition_check = parity_check_from_columns(repetition_columns, q)
    repetition_decoder = prepare_radius_four_syndrome_list_decoder(
        repetition_columns, q, parity_check=repetition_check
    )
    cosets, maximum_list_size, histogram, error_count = (
        verify_all_radius_four_syndrome_cosets(
            repetition_columns,
            repetition_decoder,
            q,
            parity_check=repetition_check,
        )
    )
    assert cosets == q**5
    assert error_count == sum(
        len(tuple(combinations(range(6), weight))) * (q - 1) ** weight
        for weight in range(5)
    )
    print(
        f"repetition [6,1,6]_{q}: syndrome-cosets={cosets}, "
        f"radius-four errors={error_count}, max-list={maximum_list_size}, "
        f"histogram={dict(sorted(histogram.items()))}"
    )

    q, k, lam = 17, 3, 1
    evaluation_set = tuple(range(7))
    source = theorem22_case_i_points(evaluation_set, k, lam, q)
    parity_check = case_i_structured_parity_check(evaluation_set, k, lam, q)
    assert k3_case_i_is_mds_by_triple_sums(evaluation_set, lam, q)
    assert is_arc(source, q)
    decoder = prepare_radius_four_syndrome_list_decoder(
        source, q, parity_check=parity_check
    )
    zero_codeword = tuple(0 for _ in source)
    support_tests = 0
    for weight in range(5):
        for support in combinations(range(len(source)), weight):
            received = [0] * len(source)
            for position, coordinate in enumerate(support):
                received[coordinate] = position + 1
            assert zero_codeword in decoder(received)
            support_tests += 1
    assert support_tests == sum(
        len(tuple(combinations(range(len(source)), weight))) for weight in range(5)
    )
    print(
        f"distance-six non-RL q={q}, k={k}, S={evaluation_set}, lambda={lam}: "
        f"radius-four support-tests={support_tests}, redundancy=5"
    )

    normalized_mds_sets = tuple(
        S
        for S in combinations(range(q), 7)
        if k3_case_i_is_mds_by_triple_sums(S, lam, q)
    )
    assert len(normalized_mds_sets) == 21
    projective_syndromes = projective_points_prime_field(5, q)
    hole_histogram: dict[int, int] = {}
    for S in normalized_mds_sets:
        check = case_i_structured_parity_check(S, k, lam, q)
        check_columns = tuple(
            tuple(check[row][coordinate] for row in range(5))
            for coordinate in range(len(check[0]))
        )
        geometric_hyperplanes = codimension_five_check_hyperplane_keys(
            check_columns, q
        )
        factorized_hyperplanes = distance_six_completion_hyperplane_keys(
            S, lam, q
        )
        assert len(set(geometric_hyperplanes)) == 126
        assert set(geometric_hyperplanes) == set(factorized_hyperplanes)
        holes = codimension_five_projective_hole_count(
            check, q, projective_points=projective_syndromes
        )
        hole_histogram[holes] = hole_histogram.get(holes, 0) + 1
    print(
        f"q={q} distance-six normalized MDS sets={len(normalized_mds_sets)}, "
        f"radius-four projective-hole histogram={dict(sorted(hole_histogram.items()))}"
    )


def run_distance_six_arrangement_scan() -> None:
    """Classify first-field holes by the first arrangement coefficients."""
    q, k, lam = 17, 3, 1
    normalized_mds_sets = tuple(
        S
        for S in combinations(range(q), 7)
        if k3_case_i_is_mds_by_triple_sums(S, lam, q)
    )
    projective_syndromes = projective_points_prime_field(5, q)
    coefficient_histogram: dict[int, dict[tuple[int, int], int]] = {}
    for S in normalized_mds_sets:
        normals = distance_six_completion_hyperplane_keys(S, lam, q)
        b2, b3, rank_two_histogram, rank_three_flat_count = (
            completion_arrangement_first_coefficients(normals, q)
        )
        assert set(rank_two_histogram) == {2, 3, 6}
        assert rank_two_histogram[6] == 84
        assert b2 == 7035 - rank_two_histogram[3]
        check = case_i_structured_parity_check(S, k, lam, q)
        holes = codimension_five_projective_hole_count(
            check, q, projective_points=projective_syndromes
        )
        by_coefficients = coefficient_histogram.setdefault(holes, {})
        by_coefficients[b2, b3] = by_coefficients.get((b2, b3), 0) + 1
        print(
            f"S={S}, holes={holes}, b2={b2}, b3={b3}, "
            f"rank2={rank_two_histogram}, rank3-flats={rank_three_flat_count}"
        )

    expected = {
        0: {
            (6746, 180747): 1,
            (6750, 181098): 2,
            (6750, 181103): 2,
            (6768, 182299): 2,
        },
        1: {
            (6737, 180262): 2,
            (6738, 180285): 2,
            (6741, 180514): 2,
            (6744, 180682): 2,
        },
        2: {
            (6744, 180735): 1,
            (6745, 180745): 2,
            (6754, 181339): 2,
        },
        4: {(6753, 181289): 1},
    }
    assert coefficient_histogram == expected
    assert len(
        {
            coefficients: holes
            for holes, entries in coefficient_histogram.items()
            for coefficients in entries
        }
    ) == 12
    print(f"distance-six arrangement classifier={coefficient_histogram}")


def run_boundary_q13_coset_scan() -> None:
    """Classify q=13 boundary cosets and verify their additive signatures."""
    q, k, evaluation_size = 13, 3, 6
    scaling_orbits: dict[
        tuple[tuple[int, ...], int], list[tuple[tuple[int, ...], int]]
    ] = {}
    for S in combinations(range(q), evaluation_size):
        for lam in range(1, q):
            arithmetic_mds = k3_boundary_case_i_is_mds(S, lam, q)
            geometric_mds = is_arc(theorem22_case_i_points(S, k, lam, q), q)
            assert arithmetic_mds == geometric_mds
            if not arithmetic_mds:
                continue
            assert sum(S) * inv_mod(lam, q) % q in (0, 1, 2)
            key = canonical_case_i_scaling_pair(S, lam, q)
            scaling_orbits.setdefault(key, []).append((S, lam))

    assert sum(map(len, scaling_orbits.values())) == 96
    assert len(scaling_orbits) == 8
    assert {len(orbit) for orbit in scaling_orbits.values()} == {q - 1}

    distribution_a = (
        (1, 853), (2, 3708), (3, 9360), (4, 9480),
        (5, 4404), (6, 684), (7, 72),
    )
    distribution_b = (
        (0, 96), (1, 481), (2, 4176), (3, 9120), (4, 9744),
        (5, 3936), (6, 1008),
    )
    distribution_c = (
        (0, 96), (1, 637), (2, 3720), (3, 9468), (4, 9804),
        (5, 3888), (6, 876), (7, 72),
    )
    distribution_d = (
        (0, 48), (1, 817), (2, 3648), (3, 9072), (4, 10416),
        (5, 3456), (6, 1104),
    )
    expected_by_signature = {
        (0, 10, 3): ("A", distribution_a),
        (1, 10, 3): ("B", distribution_b),
        (1, 11, 2): ("C", distribution_c),
        (1, 11, 3): ("D", distribution_d),
    }
    expected_normalized_types = {
        (0, 1, 2, 3, 4, 5): "A",
        (0, 1, 2, 7, 8, 9): "B",
        (0, 1, 4, 5, 6, 11): "C",
        (1, 2, 3, 4, 5, 12): "C",
        (2, 3, 6, 7, 10, 12): "D",
        (3, 4, 5, 8, 9, 11): "C",
        (4, 5, 6, 7, 8, 9): "A",
        (4, 5, 6, 7, 8, 10): "C",
    }

    distributions: dict[tuple[tuple[int, int], ...], int] = {}
    signature_pair_counts: dict[tuple[int, int, int], int] = {}
    normalized_types: dict[tuple[int, ...], str] = {}
    expected_error_count = sum(
        len(tuple(combinations(range(evaluation_size + 2), weight)))
        * (q - 1) ** weight
        for weight in range(4)
    )
    for representative, orbit in sorted(scaling_orbits.items()):
        S, lam = representative
        signature = boundary_additive_signature(S, lam, q)
        type_label, expected_distribution = expected_by_signature[signature]
        inverse_lambda = inv_mod(lam, q)
        normalized_set = tuple(sorted(a * inverse_lambda % q for a in S))
        normalized_types[normalized_set] = type_label
        parity_check = boundary_case_i_structured_parity_check(S, k, lam, q)
        cosets, errors, maximum, histogram = radius_three_syndrome_histogram(
            parity_check, q
        )
        assert cosets == q**4 and errors == expected_error_count
        distribution = tuple(histogram.items())
        assert distribution == expected_distribution
        distributions[distribution] = distributions.get(distribution, 0) + len(orbit)
        signature_pair_counts[signature] = (
            signature_pair_counts.get(signature, 0) + len(orbit)
        )
        print(
            f"normalized S={normalized_set}, type={type_label}, "
            f"signature={signature}, orbit={len(orbit)}, max-list={maximum}, "
            f"histogram={histogram}"
        )

    expected_distributions = {
        distribution_a: 24,
        distribution_b: 12,
        distribution_c: 48,
        distribution_d: 12,
    }
    assert distributions == expected_distributions
    assert signature_pair_counts == {
        (0, 10, 3): 24,
        (1, 10, 3): 12,
        (1, 11, 2): 48,
        (1, 11, 3): 12,
    }
    assert normalized_types == expected_normalized_types
    print(
        "q=13 boundary summary: MDS pairs=96, scaling orbits=8, "
        f"additive/distribution types={len(distributions)}"
    )


def run_boundary_q17_hole_scan() -> None:
    """Exhaust all normalized q=17 boundary MDS pairs and count syndrome holes."""
    q, k = 17, 3
    projective_points = projective_points_prime_field(4, q)
    hole_histogram: dict[int, int] = {}
    eta_counts = {0: 0, 1: 0}
    signature_holes: dict[tuple[int, int, int], set[int]] = {}
    normalized_mds_count = 0

    for S in combinations(range(q), 6):
        if not k3_boundary_case_i_is_mds(S, 1, q):
            continue
        normalized_mds_count += 1
        signature = boundary_additive_signature(S, 1, q)
        eta_counts[signature[0]] += 1
        parity_check = boundary_case_i_structured_parity_check(S, k, 1, q)
        holes = codimension_four_projective_hole_count(
            parity_check, q, projective_points=projective_points
        )
        hole_histogram[holes] = hole_histogram.get(holes, 0) + 1
        signature_holes.setdefault(signature, set()).add(holes)

    expected_hole_histogram = {
        22: 1, 24: 1, 29: 2, 30: 4, 31: 8, 32: 12, 33: 4,
        34: 8, 35: 12, 36: 9, 37: 18, 38: 21, 39: 14, 40: 29,
        41: 18, 42: 32, 43: 14, 44: 18, 45: 6, 46: 18, 47: 18,
        48: 18, 49: 10, 50: 14, 51: 14, 52: 18, 53: 16, 54: 4,
        57: 8,
    }
    assert normalized_mds_count == 369
    assert eta_counts == {0: 252, 1: 117}
    assert dict(sorted(hole_histogram.items())) == expected_hole_histogram
    assert min(hole_histogram) == 22 and max(hole_histogram) == 57
    assert len(signature_holes) == 19
    assert len(signature_holes[(0, 10, 3)]) > 1

    sample = (0, 1, 2, 3, 4, 5)
    sample_check = boundary_case_i_structured_parity_check(sample, k, 1, q)
    cosets, errors, maximum, histogram = radius_three_syndrome_histogram(
        sample_check, q
    )
    assert cosets == q**4 and errors == 236673 and maximum == 7
    assert histogram == {
        0: 704, 1: 7233, 2: 23696, 3: 30704,
        4: 16528, 5: 4176, 6: 416, 7: 64,
    }
    assert histogram[0] == (q - 1) * 44

    print(
        "q=17 boundary summary: normalized MDS scaling orbits=369, "
        "all covering radii=4, additive signatures=19, "
        f"projective-hole range={min(hole_histogram)}..{max(hole_histogram)}"
    )
    print(f"projective-hole histogram={dict(sorted(hole_histogram.items()))}")
    print(
        f"sample S={sample}, lambda=1: max-list={maximum}, "
        f"syndrome histogram={histogram}"
    )


def run_minimal_source_witnesses() -> None:
    """Recover source witnesses and test unique/list decoding on RL and non-RL cases."""
    q, k = 11, 3
    cases = (
        ((1, 2, 5, 6, 9), 3),
        ((3, 4, 5, 6, 7), 10),
        ((3, 4, 5, 6, 7), 8),
        ((3, 4, 5, 6, 7), 9),
    )
    expected_first_witness = {
        3: ((0, 3, 4, 7, 10), 8),
        10: ((0, 3, 4, 5, 8), 10),
        8: ((0, 4, 6, 9, 10), 7),
        9: ((0, 3, 4, 7, 10), 8),
    }
    for S, lam in cases:
        tangent_pairs = minimal_global_rl_tangent_pairs(S, k, lam, q)
        assert tangent_pairs
        first = tangent_pairs[0]
        assert (first["target_S"], first["delta"]) == expected_first_witness[lam]
        probe = tuple(range(k + 4))
        assert transport_unique_decoder(
            probe, first["monomial_map"], lambda word: word, q
        ) == probe
        source = theorem22_case_i_points(S, k, lam, q)
        target = roth_lempel_points(first["target_S"], k, first["delta"], q)
        target_decoder = prepare_single_error_syndrome_decoder(target, q)
        target_list_decoder = prepare_radius_two_syndrome_list_decoder(target, q)
        source_codeword = encode_from_columns((1, 2, 3, 4), source, q)
        assert transport_unique_decoder(
            source_codeword, first["monomial_map"], target_decoder, q
        ) == source_codeword
        for error_coordinate in range(len(source_codeword)):
            for error_value in range(1, q):
                received = list(source_codeword)
                received[error_coordinate] = (
                    received[error_coordinate] + error_value
                ) % q
                decoded = transport_unique_decoder(
                    received, first["monomial_map"], target_decoder, q
                )
                assert decoded == source_codeword
        radius_two_tests, maximum_list_size = verify_all_radius_two_error_patterns(
            source,
            source_codeword,
            lambda word: transport_list_decoder(
                word, first["monomial_map"], target_list_decoder, q
            ),
            q,
        )
        print(
            f"S={S}, lambda={lam}: S'={first['target_S']}, delta={first['delta']}, "
            f"excluded={first['excluded_label']}, tangent={first['tangency_point']}, "
            f"unique-decoder tests={1 + len(source_codeword) * (q - 1)}, "
            f"radius-two-list tests={radius_two_tests}, max-list={maximum_list_size}"
        )
        print("projectivity=", first["projectivity"])
        print("monomial_map=", first["monomial_map"])

    q, k = 13, 3
    non_rl_S, non_rl_lambda = (0, 1, 2, 3, 4), 2
    non_rl_source = theorem22_case_i_points(non_rl_S, k, non_rl_lambda, q)
    assert is_arc(non_rl_source, q)
    assert not minimal_global_rl_tangent_pairs(non_rl_S, k, non_rl_lambda, q)
    non_rl_codeword = encode_from_columns((1, 2, 3, 4), non_rl_source, q)
    non_rl_tests, non_rl_maximum_list_size = verify_all_radius_two_error_patterns(
        non_rl_source,
        non_rl_codeword,
        prepare_radius_two_syndrome_list_decoder(non_rl_source, q),
        q,
    )
    assert non_rl_tests == 3109
    non_rl_cosets, non_rl_coset_maximum_list_size = verify_all_radius_two_syndrome_cosets(
        non_rl_source,
        prepare_radius_two_syndrome_list_decoder(non_rl_source, q),
        q,
    )
    assert non_rl_cosets == q**3
    assert non_rl_coset_maximum_list_size == non_rl_maximum_list_size
    print(
        f"global non-RL S={non_rl_S}, lambda={non_rl_lambda}: "
        f"radius-two-list tests={non_rl_tests}, syndrome-cosets={non_rl_cosets}, "
        f"max-list={non_rl_maximum_list_size}"
    )


def run_global_mds_minimal_scan(q: int) -> None:
    """Exhaust the global MDS minimal-length criterion for k=3."""
    if q not in (11, 13, 17):
        raise ValueError("the verified regression counts are available for q=11,13,17")
    k, n = 3, 5
    mds_count = positive_count = tangent_pair_count = 0
    for S in combinations(range(q), n):
        for lam in range(1, q):
            points = theorem22_case_i_points(S, k, lam, q)
            geometric_mds = is_arc(points, q)
            assert geometric_mds == k3_case_i_is_mds(S, lam, q)
            if not geometric_mds:
                continue
            mds_count += 1
            tangent_pairs = minimal_global_rl_tangent_pairs(S, k, lam, q)
            tangent_pair_count += len(tangent_pairs)
            positive_count += bool(tangent_pairs)
    expected = {
        11: (210, 210, 360),
        13: (1512, 540, 720),
        17: (20592, 3600, 4800),
    }
    assert (mds_count, positive_count, tangent_pair_count) == expected[q]
    print(
        f"q={q}, k={k}, n={n}: MDS pairs={mds_count}, "
        f"global RL pairs={positive_count}, tangent_pairs={tangent_pair_count}"
    )


def run_closed_mds_minimal_scan(q: int) -> None:
    """Fast k=3 scan using only the closed MDS and RL polynomials."""
    if q < 5 or any(q % divisor == 0 for divisor in range(2, int(q**0.5) + 1)):
        raise ValueError("the self-contained implementation requires a prime field")
    mds_count = positive_count = tangent_pair_count = involution_fixed_positive = 0
    role_counts = {"evaluation": 0, "P_0": 0, "P_lambda": 0}
    rooted_translation_classes: set[tuple[int, ...]] = set()
    for S in combinations(range(q), 5):
        triple_sums = {
            (S[i] + S[j] + S[h]) % q
            for i in range(5)
            for j in range(i + 1, 5)
            for h in range(j + 1, 5)
        }
        if 0 in triple_sums:
            continue
        mds_count += q - 1 - len(triple_sums)
        total_1 = sum(S) % q
        for tangent_value in S:
            lam = (total_1 + tangent_value) % q
            if lam == 0 or lam in triple_sums:
                continue
            labels = k3_minimal_closed_excluded_labels(S, lam, q)
            candidates, cubic = k3_unified_candidate_polynomial(S, lam, q)
            assert len(candidates) == 6
            assert len({value for _, value in candidates}) == 6
            unified_labels = tuple(
                label for label, value in candidates if evaluate_cubic(cubic, value, q) == 0
            )
            assert set(unified_labels) == set(labels)
            assert len(labels) <= 2
            candidate_values = tuple(value for _, value in candidates)
            if unified_labels:
                assert not has_balanced_triple_partition(candidate_values, q)
            if not labels:
                continue
            positive_count += 1
            tangent_pair_count += len(labels)
            if q != 3:
                shift = lam * inv_mod(3, q) % q
                if tuple(sorted((shift - a) % q for a in S)) == S:
                    involution_fixed_positive += 1
            for label in labels:
                role = "evaluation" if label.startswith("evaluation_") else label
                role_counts[role] += 1
                root = next(value for candidate_label, value in candidates if candidate_label == label)
                rooted_translation_classes.add(
                    tuple(sorted((value - root) % q for value in candidate_values))
                )
            if q != 3:
                shift = lam * inv_mod(3, q) % q
                transformed_S = tuple(sorted((shift - a) % q for a in S))
                transformed_labels = set(
                    k3_minimal_closed_excluded_labels(transformed_S, lam, q)
                )
                expected_labels = {
                    "P_lambda"
                    if label == "P_0"
                    else "P_0"
                    if label == "P_lambda"
                    else f"evaluation_{(shift - int(label.removeprefix('evaluation_'))) % q}"
                    for label in labels
                }
                assert transformed_labels == expected_labels
    if q not in (2, 3):
        assert role_counts["evaluation"] == 4 * role_counts["P_0"]
        assert role_counts["evaluation"] == 4 * role_counts["P_lambda"]
        assert tangent_pair_count == 30 * len(rooted_translation_classes)
    normalized_double_roots = None
    subgroup_orbits = None
    if q > 5:
        assert len(rooted_translation_classes) == k3_closed_rooted_translation_count(q)
        normalized_double_roots = k3_normalized_double_root_count(q)
        assert normalized_double_roots == k3_closed_normalized_double_root_count(q)
        assert tangent_pair_count - positive_count == 15 * (q - 1) * normalized_double_roots
        assert positive_count == k3_closed_global_rl_pair_count(q)
        assert involution_fixed_positive == k3_closed_involution_fixed_rl_pair_count(q)
        subgroup_orbits = (positive_count + involution_fixed_positive) // (2 * (q - 1))
        assert subgroup_orbits == k3_closed_monomial_subgroup_orbit_count(q)
    known = {
        11: (210, 210, 360),
        13: (1512, 540, 720),
        17: (20592, 3600, 4800),
    }
    if q in known:
        assert (mds_count, positive_count, tangent_pair_count) == known[q]
    print(
        f"q={q}, k=3, n=5: MDS pairs={mds_count}, global RL pairs={positive_count}, "
        f"tangent_pairs={tangent_pair_count}, roles={role_counts}, "
        f"rooted translation classes={len(rooted_translation_classes)}, "
        f"normalized double-root sets={normalized_double_roots}, "
        f"involution-fixed RL pairs={involution_fixed_positive}, "
        f"subgroup orbits={subgroup_orbits}"
    )


def run_double_root_scan(q: int) -> None:
    """Evaluate the normalized overlap term and the resulting RL-pair count."""
    normalized_count = k3_normalized_double_root_count(q)
    closed_count = k3_closed_normalized_double_root_count(q)
    assert normalized_count == closed_count
    trace = k3_cm_minus8_trace(q)
    trace_square = k3_cm_minus8_trace_square(q)
    assert trace**2 == trace_square
    double_pairs = 15 * (q - 1) * normalized_count
    tangent_pair_count = 30 * k3_closed_rooted_translation_count(q)
    print(
        f"q={q}: normalized double-root sets={normalized_count}, "
        f"closed formula={closed_count}, CM trace={trace}, "
        f"CM trace square={trace_square}, "
        f"double-tangent-pair presentations={double_pairs}, "
        f"tangent_pairs={tangent_pair_count}, "
        f"global RL pairs={tangent_pair_count - double_pairs}"
    )


def run_rooted_six_set_scan(q: int) -> None:
    """Independently verify the rooted-translation-orbit tangent-pair theorem."""
    if q <= 3 or any(q % divisor == 0 for divisor in range(2, int(q**0.5) + 1)):
        raise ValueError("the rooted six-set scan requires a prime field of characteristic >3")
    inv2, inv3 = inv_mod(2, q), inv_mod(3, q)
    rooted_six_sets = 0
    rooted_classes: set[tuple[int, ...]] = set()
    tangent_pairs: set[tuple[tuple[int, ...], int, str]] = set()
    for C in combinations(range(q), 6):
        cubic = symmetric_six_set_cubic(C, q)
        sigma = sum(C) % q
        for root in C:
            if evaluate_cubic(cubic, root, q) != 0:
                continue
            rooted_six_sets += 1
            assert not has_balanced_triple_partition(C, q)
            rooted_classes.add(tuple(sorted((value - root) % q for value in C)))
            for p, z in permutations(C, 2):
                shift = (p + 2 * z - sigma) * inv3 % q
                u, v = (p + shift) % q, (z + shift) % q
                A = {(value + shift) % q for value in C if value not in (p, z)}
                m = (-u * inv2) % q
                assert len(A) == 4 and m not in A
                S = tuple(sorted(A | {m}))
                lam = (v - u) % q
                assert k3_case_i_is_mds(S, lam, q)
                candidates, source_cubic = k3_unified_candidate_polynomial(S, lam, q)
                shifted_root = (root + shift) % q
                label = next(label for label, value in candidates if value == shifted_root)
                assert evaluate_cubic(source_cubic, shifted_root, q) == 0
                tangent_pairs.add((S, lam, label))
    role_counts = {"evaluation": 0, "P_0": 0, "P_lambda": 0}
    for _, _, label in tangent_pairs:
        role = "evaluation" if label.startswith("evaluation_") else label
        role_counts[role] += 1
    assert len(tangent_pairs) == 30 * len(rooted_classes)
    assert role_counts["evaluation"] == 4 * role_counts["P_0"]
    assert role_counts["evaluation"] == 4 * role_counts["P_lambda"]
    assert len(rooted_classes) == k3_closed_rooted_translation_count(q)
    print(
        f"q={q}: rooted six-sets={rooted_six_sets}, "
        f"rooted translation classes={len(rooted_classes)}, "
        f"tangent_pairs={len(tangent_pairs)}, roles={role_counts}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--nonzero-branch-example",
        action="store_true",
        help="verify the MDS example at the degree-gap threshold",
    )
    parser.add_argument(
        "--boundary-scan",
        action="store_true",
        help="run the q=11 fixed-S regression sample on six evaluation sets",
    )
    parser.add_argument(
        "--global-boundary-scan",
        action="store_true",
        help="exhaust the global MDS criterion for k=3, n=6",
    )
    parser.add_argument(
        "--global-boundary-q",
        type=int,
        choices=(13, 17),
        default=13,
        help="prime field for --global-boundary-scan (default: 13)",
    )
    parser.add_argument(
        "--boundary-radius-three-decoder",
        action="store_true",
        help="verify the structured boundary dual and all radius-three syndrome cosets",
    )
    parser.add_argument(
        "--distance-six-radius-four-decoder",
        action="store_true",
        help="verify exact radius-four decoding at redundancy five",
    )
    parser.add_argument(
        "--distance-six-arrangement-scan",
        action="store_true",
        help="classify q=17 distance-six holes by arrangement coefficients",
    )
    parser.add_argument(
        "--boundary-q13-coset-scan",
        action="store_true",
        help="classify radius-three coset distributions on all q=13 boundary MDS pairs",
    )
    parser.add_argument(
        "--boundary-q17-hole-scan",
        action="store_true",
        help="exhaust projective syndrome holes on all normalized q=17 boundary MDS pairs",
    )
    parser.add_argument(
        "--minimal-source-witnesses",
        action="store_true",
        help="recover unrestricted RL witnesses for the minimal source-paper cases",
    )
    parser.add_argument(
        "--global-minimal-scan",
        action="store_true",
        help="exhaust the global MDS criterion for k=3, n=5",
    )
    parser.add_argument(
        "--global-minimal-q",
        type=int,
        choices=(11, 13, 17),
        default=11,
        help="prime field for --global-minimal-scan (default: 11)",
    )
    parser.add_argument(
        "--closed-minimal-scan",
        action="store_true",
        help="run the fast closed-polynomial k=3 minimal scan",
    )
    parser.add_argument(
        "--closed-minimal-q",
        type=int,
        default=19,
        help="prime field for --closed-minimal-scan (default: 19)",
    )
    parser.add_argument(
        "--rooted-six-scan",
        action="store_true",
        help="independently verify the rooted six-set translation-orbit theorem",
    )
    parser.add_argument(
        "--rooted-six-q",
        type=int,
        default=13,
        help="prime field for --rooted-six-scan (default: 13)",
    )
    parser.add_argument(
        "--double-root-scan",
        action="store_true",
        help="evaluate the normalized double-tangent-pair overlap term",
    )
    parser.add_argument(
        "--double-root-q",
        type=int,
        default=43,
        help="prime field for --double-root-scan (default: 43)",
    )
    args = parser.parse_args()
    if args.nonzero_branch_example:
        run_nonzero_branch_threshold_example()
    elif args.boundary_q17_hole_scan:
        run_boundary_q17_hole_scan()
    elif args.boundary_q13_coset_scan:
        run_boundary_q13_coset_scan()
    elif args.distance_six_arrangement_scan:
        run_distance_six_arrangement_scan()
    elif args.distance_six_radius_four_decoder:
        run_distance_six_radius_four_decoder()
    elif args.boundary_radius_three_decoder:
        run_boundary_radius_three_decoder()
    elif args.minimal_source_witnesses:
        run_minimal_source_witnesses()
    elif args.double_root_scan:
        run_double_root_scan(args.double_root_q)
    elif args.rooted_six_scan:
        run_rooted_six_set_scan(args.rooted_six_q)
    elif args.closed_minimal_scan:
        run_closed_mds_minimal_scan(args.closed_minimal_q)
    elif args.global_minimal_scan:
        run_global_mds_minimal_scan(args.global_minimal_q)
    elif args.global_boundary_scan:
        run_global_mds_boundary_scan(args.global_boundary_q)
    elif args.boundary_scan:
        run_q11_boundary_scan()
    else:
        run_nonzero_branch_threshold_example()
        run_closed_mds_minimal_scan(13)
        run_rooted_six_set_scan(13)
        run_double_root_scan(23)
