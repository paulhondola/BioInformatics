"""Bonus: matplotlib drawing of the rooted NJ tree, leaf labels coloured by expected group."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from matplotlib.patches import Patch

from phylogenetics.dataset import SPECIES

GROUP_COLORS: dict[str, str] = {
    "Primates": "#c0392b",
    "Rodents": "#d68910",
    "Laurasiatheria": "#7d3c98",
    "Bird": "#1e8449",
    "Amphibian": "#2471a3",
    "Fish": "#566573",
}


def draw_tree(tree: Tree, out: Path) -> Path:
    """Save a branch-length-scaled phylogram of ``tree`` with a group legend to ``out``."""
    color_of = {s.label: GROUP_COLORS[s.group] for s in SPECIES}
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=lambda c: c.name.replace("_", " ") if c.is_terminal() else None,
        label_colors=lambda name: color_of.get(name.replace(" ", "_"), "black"),
    )
    ax.set_title("TP53 CDS Neighbor-Joining tree (Jukes-Cantor), fish-rooted")
    ax.set_xlabel("Substitutions per site")
    ax.set_ylabel("")
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.legend(
        handles=[Patch(color=c, label=g) for g, c in GROUP_COLORS.items()],
        loc="upper left",
        frameon=False,
        fontsize=8,
    )
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out
