"""
Enumeration of pattern-based state-conserving cellular automata of type 0
with the four-cell neighbourhood N = (-1, 0, 1, 2).

A pattern-based CA F_T is defined by an admissible set T of elementary
transformations, T = (conditional swaps) U (bypasses):

    conditional swap  <p,q>_c : pattern 'p q c'  ->  'q p c'           (p != q)
    bypass            [p,q,b] : pattern 'p q b'  ->  'b p q'
                                pattern '~p q b' ->  '~p b q'           (p != b, q != b)

The script enumerates three disjoint families of admissible sets:
    1. conditional swaps only (including the empty set, i.e. the identity),
    2. bypasses only,
    3. mixed sets (at least one bypass and at least one conditional swap),
builds the lookup table of every resulting local rule, optionally verifies
state conservation with the Boccara-Fuks criterion, and writes the results
to text files.

Usage:
    python code.py -k 3 -o results --verify
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import chain, combinations, product
from math import comb
from pathlib import Path
from typing import Dict, FrozenSet, Iterable, Iterator, List, NamedTuple, Sequence, Set, Tuple, TypeVar

State = int
Neighbourhood = Tuple[State, State, State, State]
LocalRule = Dict[Neighbourhood, State]
NEIGHBOURHOOD_SIZE = 4


# ---------------------------------------------------------------------------
# Elementary transformations
# ---------------------------------------------------------------------------

class ConditionalSwap(NamedTuple):
    """Conditional swap <p,q>_c: pattern 'p q c' becomes 'q p c'."""
    p: State
    q: State
    c: State

    def __str__(self) -> str:
        return f"<{self.p},{self.q}>_{self.c}"


class Bypass(NamedTuple):
    """Bypass [p,q,b]: 'p q b' becomes 'b p q' and '~p q b' becomes '~p b q'.
    p, q are the bypassed states, b is the bypassing state."""
    p: State
    q: State
    b: State

    def __str__(self) -> str:
        return f"[{self.p},{self.q},{self.b}]"


SwapSet = Tuple[ConditionalSwap, ...]
BypassSet = Tuple[Bypass, ...]


class TransformationSet(NamedTuple):
    """Admissible set T of elementary transformations defining F_T."""
    bypasses: BypassSet
    swaps: SwapSet

    def __str__(self) -> str:
        """Bypasses followed by conditional swaps, separated by spaces (empty for T = {})."""
        return " ".join(map(str, (*self.bypasses, *self.swaps)))


def all_conditional_swaps(k: int) -> Tuple[ConditionalSwap, ...]:
    """The set S_k of all conditional swaps, |S_k| = k^2 (k-1)."""
    return tuple(ConditionalSwap(p, q, c)
                 for p, q, c in product(range(k), repeat=3) if p != q)


def all_bypasses(k: int) -> Tuple[Bypass, ...]:
    """The set B_k of all bypasses, |B_k| = k (k-1)^2."""
    return tuple(Bypass(p, q, b)
                 for p, q, b in product(range(k), repeat=3) if p != b and q != b)


# ---------------------------------------------------------------------------
# Conflicts
# ---------------------------------------------------------------------------

def swaps_in_conflict(s: ConditionalSwap, t: ConditionalSwap) -> bool:
    """<a,b>_c is in conflict with <b,c>_* and with <~a,a>_b.
    The second case is the first one read in the opposite direction."""
    def blocks(x: ConditionalSwap, y: ConditionalSwap) -> bool:
        return y.p == x.q and y.q == x.c          # y = <b,c>_* for x = <a,b>_c
    return blocks(s, t) or blocks(t, s)


def bypasses_in_conflict(u: Bypass, v: Bypass) -> bool:
    """Two bypasses conflict if a state is bypassing for one and bypassed for the other."""
    return u.b in (v.p, v.q) or v.b in (u.p, u.q)


def bypass_and_swap_in_conflict(bypass: Bypass, swap: ConditionalSwap) -> bool:
    """Conflicts of [p,q,b] with conditional swaps."""
    p, q, b = bypass
    return (
        (swap.q == p and swap.c == q)       # <~p,p>_q
        or (swap.q == q and swap.c == b)    # <*,q>_b, in particular <p,q>_b
        or (swap.p == q and swap.q == b)    # <q,b>_*
        or swap.p == b                      # <b,~b>_*
    )


def conflict_graph(swaps: Sequence[ConditionalSwap]) -> Dict[ConditionalSwap, Set[ConditionalSwap]]:
    """Conflict graph Gamma_{S'} as an adjacency map."""
    graph: Dict[ConditionalSwap, Set[ConditionalSwap]] = {s: set() for s in swaps}
    for s, t in combinations(swaps, 2):
        if swaps_in_conflict(s, t):
            graph[s].add(t)
            graph[t].add(s)
    return graph


X = TypeVar("X")


def independent_sets(graph: Dict[X, Set[X]], include_empty: bool = True) -> Iterator[Tuple[X, ...]]:
    """All independent sets of a graph (backtracking over the vertex order).
    With include_empty=False the empty set is never constructed: every branch
    starts by fixing the first chosen vertex."""
    vertices = list(graph)
    index = {v: i for i, v in enumerate(vertices)}
    neighbours = [{index[u] for u in graph[v]} for v in vertices]
    n = len(vertices)

    def extend(i: int, chosen: List[int], blocked: Set[int]) -> Iterator[Tuple[X, ...]]:
        if i == n:
            yield tuple(vertices[j] for j in chosen)
            return
        yield from extend(i + 1, chosen, blocked)
        if i not in blocked:
            chosen.append(i)
            yield from extend(i + 1, chosen, blocked | neighbours[i])
            chosen.pop()

    if include_empty:
        yield ()
    for first in range(n):
        yield from extend(first + 1, [first], set(neighbours[first]))


# ---------------------------------------------------------------------------
# Admissible sets
# ---------------------------------------------------------------------------

def admissible_swap_sets(k: int, include_empty: bool = True) -> Iterator[SwapSet]:
    """Admissible subsets of S_k = independent sets of Gamma_{S_k}."""
    yield from independent_sets(conflict_graph(all_conditional_swaps(k)), include_empty)


def nonempty_subsets(items: Sequence[X]) -> Iterator[Tuple[X, ...]]:
    return chain.from_iterable(combinations(items, r) for r in range(1, len(items) + 1))


def admissible_bypass_sets(k: int) -> Iterator[BypassSet]:
    """Non-empty admissible subsets of B_k.

    Such a set is determined by its set of bypassing states B (non-empty, proper
    subset of S) and, for every b in B, a non-empty set of pairs (p, q) of
    bypassed states taken from S \\ B."""
    states = range(k)
    for size in range(1, k):
        for bypassing in combinations(states, size):
            bypassed = [s for s in states if s not in bypassing]
            pairs = list(product(bypassed, repeat=2))
            options = [[tuple(Bypass(p, q, b) for p, q in chosen) for chosen in nonempty_subsets(pairs)]
                       for b in bypassing]
            for choice in product(*options):
                yield tuple(chain.from_iterable(choice))


def compatible_swaps(k: int, bypasses: BypassSet) -> Tuple[ConditionalSwap, ...]:
    """Conditional swaps that are in conflict with none of the given bypasses."""
    return tuple(s for s in all_conditional_swaps(k)
                 if not any(bypass_and_swap_in_conflict(b, s) for b in bypasses))


def admissible_mixed_sets(k: int) -> Iterator[TransformationSet]:
    """Admissible sets containing at least one bypass and at least one conditional swap.
    The empty set of conditional swaps is never generated, and bypass sets which
    leave no compatible conditional swap are skipped."""
    for bypasses in admissible_bypass_sets(k):
        swaps = compatible_swaps(k, bypasses)
        if not swaps:
            continue
        for swap_set in independent_sets(conflict_graph(swaps), include_empty=False):
            yield TransformationSet(bypasses, swap_set)


def all_pattern_based(k: int) -> Iterator[Tuple[str, TransformationSet]]:
    """All pattern-based CAs of type 0, labelled by family."""
    for swaps in admissible_swap_sets(k, include_empty=True):
        yield "conditional_swaps", TransformationSet((), swaps)
    for bypasses in admissible_bypass_sets(k):
        yield "bypasses", TransformationSet(bypasses, ())
    for mixed in admissible_mixed_sets(k):
        yield "mixed", mixed


# ---------------------------------------------------------------------------
# Counting formulas (cross-check of the enumeration)
# ---------------------------------------------------------------------------

def count_swap_sets_formula(k: int) -> int:
    """sum over p in {0,1}^{k(k-1)} of prod_{i != j} (2^{k - p_{j,*}} - 1)^{p_{i,j}}."""
    pairs = [(i, j) for i in range(k) for j in range(k) if i != j]
    total = 0
    for bits in product((0, 1), repeat=len(pairs)):
        p = dict(zip(pairs, bits))
        out_degree = {j: sum(p[(j, l)] for l in range(k) if l != j) for j in range(k)}
        term = 1
        for (i, j), used in p.items():
            if used:
                term *= 2 ** (k - out_degree[j]) - 1
        total += term
    return total


def count_bypass_sets_formula(k: int) -> int:
    """sum_{i=1}^{k-1} C(k, i) (2^{(k-i)^2} - 1)^i  (i = number of bypassing states)."""
    return sum(comb(k, i) * (2 ** ((k - i) ** 2) - 1) ** i for i in range(1, k))


# ---------------------------------------------------------------------------
# Local rules and lookup tables
# ---------------------------------------------------------------------------

def local_rule(k: int, transformations: TransformationSet) -> LocalRule:
    """Local rule f of F_T for the neighbourhood (-1, 0, 1, 2);
    the updated cell is the second entry of the tuple."""
    states = range(k)
    f: LocalRule = {x: x[1] for x in product(states, repeat=NEIGHBOURHOOD_SIZE)}

    for p, q, c in transformations.swaps:              # 'p q c' -> 'q p c'
        for s in states:
            f[(s, p, q, c)] = q                        # cell in state p
            f[(p, q, c, s)] = p                        # cell in state q

    bypassed_by: Dict[Tuple[State, State], Set[State]] = defaultdict(set)
    for p, q, b in transformations.bypasses:           # 'p q b' -> 'b p q'
        bypassed_by[(q, b)].add(p)
        for s in states:
            f[(s, p, q, b)] = b                        # cell in state p
            f[(p, q, b, s)] = p                        # cell in state q
            for s2 in states:
                f[(q, b, s, s2)] = q                   # cell in state b

    for (q, b), ps in bypassed_by.items():             # '~p q b' -> '~p b q'
        for other in states:
            if other not in ps:
                for s in states:
                    f[(other, q, b, s)] = b            # cell in state q

    return f


def lookup_table(f: LocalRule) -> str:
    """Lookup table as a string of outputs, from (k-1,k-1,k-1,k-1) down to (0,0,0,0)."""
    return "".join(str(f[x]) for x in sorted(f, reverse=True))


def is_state_conserving(f: LocalRule, k: int) -> bool:
    """Boccara-Fuks criterion for every state alpha in {1, ..., k-1}."""
    m = NEIGHBOURHOOD_SIZE
    for x in product(range(k), repeat=m):
        for alpha in range(1, k):
            def delta(s: State) -> int:
                return int(s == alpha)
            rhs = delta(x[0]) + sum(
                delta(f[(0,) * (m - j) + x[1:j + 1]]) - delta(f[(0,) * (m - j) + x[:j]])
                for j in range(1, m)
            )
            if delta(f[x]) != rhs:
                return False
    return True


def reflection(f: LocalRule) -> LocalRule:
    """Reflected rule f^R(x1, x2, x3, x4) = f(x4, x3, x2, x1).
    By the reflection lemma, F in SC(tau) implies F^R in SC(1 - tau)."""
    return {x: f[x[::-1]] for x in f}


def shift_rules(k: int) -> Dict[str, Tuple[int, LocalRule]]:
    """The two state-conserving rules outside SC(0) and SC(1), with their types:
    shift-right f = x1 (type -1) and shift-left by two cells f = x4 (type 2)."""
    neighbourhoods = list(product(range(k), repeat=NEIGHBOURHOOD_SIZE))
    return {
        "shift-right": (-1, {x: x[0] for x in neighbourhoods}),
        "shift-left-2": (2, {x: x[3] for x in neighbourhoods}),
    }


def rule_type(f: LocalRule) -> int:
    """Type tau of the projection onto {0, 1}: tau + 1 of f(0,0,0,1), f(0,0,1,1), f(0,1,1,1) equal 1."""
    return f[(0, 0, 0, 1)] + f[(0, 0, 1, 1)] + f[(0, 1, 1, 1)] - 1


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

FAMILIES = ("conditional_swaps", "bypasses", "mixed")


class Checker:
    """Optional verification: state conservation, expected type, uniqueness."""

    def __init__(self, k: int, enabled: bool) -> None:
        self.k, self.enabled = k, enabled
        self.seen: Set[str] = set()

    def __call__(self, f: LocalRule, lut: str, expected_type: int, label: str) -> None:
        if not self.enabled:
            return
        if not is_state_conserving(f, self.k):
            raise AssertionError(f"not state-conserving: {label}")
        if rule_type(f) != expected_type:
            raise AssertionError(f"type {rule_type(f)} instead of {expected_type}: {label}")
        if lut in self.seen:
            raise AssertionError(f"duplicate lookup table: {label}")
        self.seen.add(lut)


def write_results(k: int, output_dir: Path, details: bool, verify: bool) -> Dict[str, int]:
    """Writes, for k states:
        k{k}_conditional_swaps.txt, k{k}_bypasses.txt, k{k}_mixed.txt  - the three families,
        k{k}_type0.txt              - all pattern-based CAs of type 0,
        k{k}_state_conserving.txt   - all state-conserving CAs (type 0, type 1 = reflections, two shifts).
    Every line is the lookup table followed by the elements of T (bypasses, then conditional
    swaps), all separated by spaces; for T = {} the line holds only the lookup table.
    For type 1 rules T is the set of the reflected rule, i.e. the rule is (F_T)^R.
    With details=False every line contains only the lookup table."""
    output_dir.mkdir(parents=True, exist_ok=True)
    names = (*FAMILIES, "type0", "state_conserving")
    files = {name: open(output_dir / f"k{k}_{name}.txt", "w", encoding="utf-8") for name in names}
    counts = {name: 0 for name in names}
    check = Checker(k, verify)

    def write(name: str, lut: str, *columns: str) -> None:
        line = " ".join(filter(None, (lut, *columns))) if details else lut
        files[name].write(line + "\n")
        counts[name] += 1

    try:
        # type 0: the three families of pattern-based CAs
        for family, transformations in all_pattern_based(k):
            f = local_rule(k, transformations)
            lut = lookup_table(f)
            check(f, lut, 0, str(transformations))
            construction = str(transformations)
            write(family, lut, construction)
            write("type0", lut, construction)
            write("state_conserving", lut, construction)

        # type 1: reflections of type 0 (second pass keeps the output grouped by type)
        for _, transformations in all_pattern_based(k):
            f = reflection(local_rule(k, transformations))
            lut = lookup_table(f)
            check(f, lut, 1, f"reflection of {transformations}")
            construction = str(transformations)
            write("state_conserving", lut, construction)

        # the two shifts
        for name, (tau, f) in shift_rules(k).items():
            lut = lookup_table(f)
            check(f, lut, tau, name)
            write("state_conserving", lut)
    finally:
        for handle in files.values():
            handle.close()
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description="Pattern-based state-conserving CAs, four-cell neighbourhood.")
    parser.add_argument("-k", "--states", type=int, default=3, help="number of states k")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--lut-only", action="store_true",
                        help="write only lookup tables, without transformation sets")
    parser.add_argument("--verify", action="store_true",
                        help="check state conservation, type and uniqueness of every lookup table")
    args = parser.parse_args()
    k = args.states

    counts = write_results(k, args.output_dir, details=not args.lut_only, verify=args.verify)

    print(f"k = {k}")
    print(f"conditional swaps only (incl. identity): {counts['conditional_swaps']}"
          f"   formula: {count_swap_sets_formula(k)}")
    print(f"bypasses only:                           {counts['bypasses']}"
          f"   formula: {count_bypass_sets_formula(k)}")
    print(f"mixed:                                   {counts['mixed']}")
    print(f"pattern-based CAs of type 0:             {counts['type0']}")
    print(f"state-conserving CAs: {counts['state_conserving']}")
    if args.verify:
        print("verified: all lookup tables are distinct, state-conserving and of the expected type")


if __name__ == "__main__":
    main()