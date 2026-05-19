import random

import numpy as np

from constants import NUCLEOTIDES


class SubstitutionModel:
    name: str
    matrix: np.ndarray

    def __init__(self, name, matrix):
        self.name = name
        self.matrix = np.asarray(matrix, dtype=float)
        self.validate_matrix()

    def validate_matrix(self) -> None:
        if len(self.matrix) != 4:
            raise ValueError("Transition matrix must have 4 rows")
        for row in self.matrix:
            if len(row) != 4:
                raise ValueError("Transition matrix must have 4 columns per row")
            row_sum = sum(row)
            if abs(row_sum - 1.0) > 1e-8:
                raise ValueError(f"Each row of the transition matrix must sum to 1.0, got {row_sum}")

    def substitute_nucleotide(
        self,
        ref_nt: str,
        rng: random.Random,
    ) -> str:
        """Sample the child nucleotide from this model's row for ``ref_nt``."""
        idx = NUCLEOTIDES.index(ref_nt)
        row = self.matrix[idx]
        return rng.choices(NUCLEOTIDES, weights=row, k=1)[0]


class JukesCantor(SubstitutionModel):
    """Equal substitution rates to every other nucleotide (diagonal zero, off-diagonal 1/3)."""

    def __init__(self):
        p = 1.0 / 3.0
        matrix = np.full((4, 4), p)
        np.fill_diagonal(matrix, 0.0)
        super().__init__(name="Jukes-Cantor", matrix=matrix)
