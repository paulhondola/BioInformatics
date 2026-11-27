"""Task 2: Neighbor-Joining tree, rooting, unrooted splits and monophyly checks."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

OUTGROUP: tuple[str, str] = ("Danio_rerio", "Oryzias_latipes")


@dataclass(frozen=True)
class NJResult:
    """Rooted NJ tree and how it was rooted (``"outgroup"`` or ``"midpoint"``)."""

    tree: Tree
    rooting: str


def splits(tree: Tree) -> set[frozenset[str]]:
    """Non-trivial bipartitions of the unrooted tree, each as the side without the smallest leaf name."""
    leaves = frozenset(t.name for t in tree.get_terminals())
    ref = min(leaves)
    out: set[frozenset[str]] = set()
    for clade in tree.find_clades():
        side = frozenset(t.name for t in clade.get_terminals())
        if ref in side:
            side = leaves - side
        if 2 <= len(side) <= len(leaves) - 2:
            out.add(side)
    return out


def is_split(tree: Tree, members: Iterable[str]) -> bool:
    """True if ``members`` (or its complement) is one side of an edge of the unrooted tree."""
    group = frozenset(members)
    leaves = frozenset(t.name for t in tree.get_terminals())
    return group in splits(tree) or (leaves - group) in splits(tree)


def build_nj_tree(
    labels: list[str], matrix: list[list[float]], outgroup: Sequence[str] = OUTGROUP
) -> NJResult:
    """NJ on a full square matrix; root on ``outgroup`` if it is a split, else at the midpoint."""
    dm = DistanceMatrix(list(labels), [row[: i + 1] for i, row in enumerate(matrix)])
    tree = DistanceTreeConstructor().nj(dm)
    if is_split(tree, outgroup):
        # Root on an ingroup leaf first so the outgroup becomes a proper clade, then put the root
        # halfway along that clade's stem. Rooting on the leaf list directly makes the outgroup's own
        # node the root, leaving a trifurcation (outgroup leaf, outgroup leaf, rest).
        ingroup = next(t for t in tree.get_terminals() if t.name not in outgroup)
        tree.root_with_outgroup(ingroup)
        stem = tree.common_ancestor(*outgroup)
        tree.root_with_outgroup(
            stem, outgroup_branch_length=(stem.branch_length or 0) / 2
        )
        rooting = "outgroup"
    else:
        print(
            f"Warning: outgroup {', '.join(outgroup)} is not a clade in the NJ tree; "
            "rooting at the midpoint instead"
        )
        tree.root_at_midpoint()
        rooting = "midpoint"
    for clade in tree.get_nonterminals():
        clade.name = None
    return NJResult(tree, rooting)


def monophyly(tree: Tree, groups: dict[str, frozenset[str]]) -> dict[str, bool]:
    """For every group with at least 2 members: is it a split of the unrooted tree?"""
    return {name: is_split(tree, m) for name, m in groups.items() if len(m) >= 2}


def write_newick(tree: Tree, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, path, "newick")


def format_task2_report(result: NJResult, checks: dict[str, bool]) -> str:
    """ASCII tree, rooting method, and one yes/no line per expected group."""
    buf = StringIO()
    Phylo.draw_ascii(result.tree, file=buf)
    rooting = (
        f"rooted on the fish outgroup ({', '.join(OUTGROUP)})"
        if result.rooting == "outgroup"
        else "midpoint-rooted (fish outgroup was not a clade)"
    )
    return "\n".join(
        [
            "=== TASK 2: Neighbor-Joining tree (Jukes-Cantor distances) ===",
            buf.getvalue().rstrip(),
            f"Tree is {rooting}.",
            "Expected groups recovered as clades (unrooted splits):",
            *(f"  {name}: {'yes' if ok else 'no'}" for name, ok in checks.items()),
        ]
    )
