import math
import random

from substitution_model import SubstitutionModel


def poisson(rng: random.Random, lam: float) -> int:
    if lam <= 0.0:
        return 0
    L = math.exp(-lam)
    k = 0
    p = 1.0
    while p > L:
        k += 1
        p *= rng.random()
    return k - 1


def molecular_clock_evolve(
    seq: str,
    branch_time: float,
    rate: float,
    model: SubstitutionModel,
    rng: random.Random,
) -> str:
    """
    Evolve a sequence along a branch of length ``branch_time``.

    Each site accumulates Poisson(rate * branch_time) substitution events; each event
    replaces the nucleotide using that row of the substitution model matrix.
    ``rate`` must use the same time units as ``branch_time`` (e.g. per site per simulation day).
    """
    if branch_time <= 0.0 or rate <= 0.0:
        return seq
    out = list(seq)
    lam = rate * branch_time
    for i, nt in enumerate(out):
        n_events = poisson(rng, lam)
        for _ in range(n_events):
            nt = model.substitute_nucleotide(ref_nt=nt, rng=rng)
        out[i] = nt
    return "".join(out)
