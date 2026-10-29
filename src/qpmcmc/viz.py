"""Figures: one accent colour plus greys."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

ACCENT = "#C2410C"
DARK, MID, LIGHT = "#222222", "#777777", "#BBBBBB"
CMAP = LinearSegmentedColormap.from_list("grey_accent", ["#F4F4F4", "#D9A38C", ACCENT])
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#E5E5E5", "grid.linewidth": 0.6})


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def waterfall(labels: list[str], start: float, steps: list[float], path: Path, ylabel: str, title: str) -> None:
    """Waterfall: first bar is the unmitigated bias, then each stage's (negative) contribution, then the residual."""
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    level = start
    ax.bar(0, start, color=DARK)
    ax.text(0, start, f"{start:.3g}", ha="center", va="bottom", fontsize=9)
    for i, d in enumerate(steps, start=1):
        ax.bar(i, d, bottom=level, color=ACCENT if d < 0 else MID)
        ax.plot([i - 1 + 0.4, i - 0.4], [level, level], color=MID, lw=0.8)
        level += d
        ax.text(i, max(level, level - d), f"{d:+.2g}", ha="center", va="bottom", fontsize=9)
    ax.bar(len(steps) + 1, level, color=LIGHT)
    ax.text(len(steps) + 1, level, f"{level:.2g}", ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, start * 1.12)
    ax.set_xticks(range(len(steps) + 2))
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="x", visible=False)
    _save(fig, path)


def line_band(x, series: dict[str, tuple[np.ndarray, np.ndarray]], path: Path, xlabel: str, ylabel: str,
              title: str, hline: tuple[float, str] | None = None, logy: bool = True) -> None:
    """Mean line with a shaded mean +/- std band for each series."""
    styles = [(DARK, "-", "o"), (ACCENT, "--", "s"), (MID, "-.", "^"), (LIGHT, ":", "D")]
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    for (name, (m, s)), (c, ls, mk) in zip(series.items(), styles):
        m, s = np.asarray(m), np.asarray(s)
        ax.plot(x, m, color=c, ls=ls, marker=mk, label=name, lw=1.8, ms=5)
        ax.fill_between(x, np.clip(m - s, 1e-12 if logy else None, None), m + s, color=c, alpha=0.18, lw=0)
    if hline is not None:
        ax.axhline(hline[0], color=MID, lw=1, ls=(0, (1, 2)))
        ax.text(x[0], hline[0] * 1.12, hline[1], color=MID, fontsize=8.5, va="bottom")
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=9)
    _save(fig, path)


def heatmap(values: np.ndarray, xlabels, ylabels, path: Path, xlabel: str, ylabel: str, title: str,
            cbar: str, fmt: str = "{:.3f}") -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    im = ax.imshow(values, cmap=CMAP, aspect="auto", origin="lower")
    ax.set_xticks(range(len(xlabels)), [str(v) for v in xlabels])
    ax.set_yticks(range(len(ylabels)), [str(v) for v in ylabels])
    thr = values.max() * 0.6
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            ax.text(j, i, fmt.format(values[i, j]), ha="center", va="center", fontsize=8.5,
                    color="white" if values[i, j] > thr else DARK)
    ax.grid(False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=cbar)
    _save(fig, path)
