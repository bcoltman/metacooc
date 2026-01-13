#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_analysis_obj(
    df: pd.DataFrame,
    out_file: str,
    q_thresh: float = None,
):
    """
    Four fixed panels (all required columns assumed present):

        Panel 1: phi vs taxon index (sorted by p_T_given_X)
        Panel 2: p_T_given_X vs taxon index
        Panel 3: phi vs log2(effect)
        Panel 4: phi vs log10(chi2)

    Required columns:
        phi, p_T_given_X, chi2,
        log2_RR_T or log2_RR_X or RR_T or RR_X, q_bh optionally
    """

    if df.empty:
        raise ValueError("DataFrame is empty — nothing to plot.")

    # Clean infinities
    df = df.replace([np.inf, -np.inf], np.nan)

    # -------- Sorting (always by specificity p_T_given_X) ----------
    df_sorted = df.sort_values("p_T_given_X", ascending=False)
    # df_sorted = df.sort_values("phi", ascending=False)
    x_idx = np.arange(len(df_sorted))

    phi_vals = df_sorted["phi"].values
    pTx_vals = df_sorted["p_T_given_X"].values
    chi_vals = df_sorted["chi2"].values

    # Colours by phi sign
    colours = np.where(phi_vals >= 0, "#1f77b4", "#d62728")

    # -------- Effect size: pick first available ----------
    if "log2_RR_T" in df_sorted.columns:
        effect = df_sorted["log2_RR_T"].values
        effect_label = "log2 RR_T"
    elif "log2_RR_X" in df_sorted.columns:
        effect = df_sorted["log2_RR_X"].values
        effect_label = "log2 RR_X"
    elif "RR_T" in df_sorted.columns:
        effect = np.log2(df_sorted["RR_T"].replace(0, np.nan).values)
        effect_label = "log2 RR_T"
    else:
        effect = np.log2(df_sorted["RR_X"].replace(0, np.nan).values)
        effect_label = "log2 RR_X"

    log_chi = np.log10(chi_vals + 1e-12)

    # -------- Set up 4-panel figure ----------
    # fig, axes = plt.subplots(
        # 4, 1,
        # figsize=(10, 16),
        # gridspec_kw={"height_ratios": [2, 2, 2, 2]}
    # )
    fig, axes = plt.subplots(
        2, 1,
        figsize=(10, 10),
        gridspec_kw={"height_ratios": [2, 2]}
    )

    # ============================================================
    # Panel 1 — phi vs index (sorted by p_T_given_X)
    # ============================================================
    ax = axes[0]
    ax.scatter(x_idx, phi_vals, c=colours, s=30, alpha=0.8)
    ax.axhline(0, color="black", lw=1)
    ax.set_ylabel("phi", fontsize=14)
    ax.set_title("Association Strength (phi)", fontsize=16)
    ax.grid(True, linestyle="--", alpha=0.6)

    if q_thresh is not None and "q_bh" in df_sorted.columns:
        sig = df_sorted["q_bh"] <= q_thresh
        if sig.any():
            ax.scatter(
                x_idx[sig.values],
                phi_vals[sig.values],
                color="gold",
                s=60,
                edgecolor="black",
                label=f"q ≤ {q_thresh}"
            )
            ax.legend()

    ax.set_xlabel("Taxa (sorted by p_T_given_X)", fontsize=14)

    # ============================================================
    # Panel 2 — p_T_given_X vs index
    # ============================================================
    ax = axes[1]
    ax.scatter(x_idx, pTx_vals, c=colours, s=30, alpha=0.8)
    ax.set_ylabel("p_T_given_X", fontsize=14)
    ax.set_title("Specificity p(T|X) Across Taxa", fontsize=16)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.set_xlabel("Taxa (sorted by p_T_given_X)", fontsize=14)

    # # ============================================================
    # # Panel 3 — phi vs log2(effect)
    # # ============================================================
    # ax = axes[2]
    # finite_mask = np.isfinite(effect) & np.isfinite(phi_vals)
    # ax.scatter(
        # effect[finite_mask],
        # phi_vals[finite_mask],
        # c=np.where(phi_vals[finite_mask] >= 0, "#1f77b4", "#d62728"),
        # s=30,
        # alpha=0.8,
    # )
    # ax.axvline(0, color="black", lw=1)
    # ax.axhline(0, color="black", lw=1)
    # ax.set_xlabel(effect_label, fontsize=14)
    # ax.set_ylabel("phi", fontsize=14)
    # ax.set_title("Association vs Effect Size", fontsize=16)
    # ax.grid(True, linestyle="--", alpha=0.6)

    # # ============================================================
    # # Panel 4 — phi vs log10(chi2)
    # # ============================================================
    # ax = axes[3]
    # finite_mask = np.isfinite(log_chi) & np.isfinite(phi_vals)
    # ax.scatter(
        # log_chi[finite_mask],
        # phi_vals[finite_mask],
        # c=np.where(phi_vals[finite_mask] >= 0, "#1f77b4", "#d62728"),
        # s=30,
        # alpha=0.8,
    # )
    # ax.axvline(0, color="black", lw=1)
    # ax.axhline(0, color="black", lw=1)
    # ax.set_xlabel("log10(chi2)", fontsize=14)
    # ax.set_ylabel("phi", fontsize=14)
    # ax.set_title("Association vs χ²", fontsize=16)
    # ax.grid(True, linestyle="--", alpha=0.6)

    # Save
    fig.tight_layout()
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[plot_analysis] Saved: {out_file}")


def plot_analysis(
    df_file: str,
    output_dir: str,
    tag: str = "",
    q_thresh: float = None,
):
    if not os.path.exists(df_file):
        raise FileNotFoundError(df_file)

    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(df_file, sep="\t")
    out_file = os.path.join(output_dir, f"{tag}phi_effect.png")

    plot_analysis_obj(df, out_file, q_thresh=q_thresh)