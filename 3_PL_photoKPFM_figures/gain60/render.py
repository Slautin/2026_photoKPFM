"""Render six individual gain-60 panels from recomputed public source tables."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

LABEL_WELLS = ("G1", "B11", "C12", "H11", "D5")
TRACE_COLORS = {"G1": "#1b9e77", "B11": "#7570b3", "C12": "#66a61e", "H11": "#a6761d", "D5": "#666666"}
HEIGHT = np.sqrt(3) / 2


def style():
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 11, 'axes.labelsize': 12, 'axes.linewidth': 1.0, 'xtick.labelsize': 10, 'ytick.labelsize': 10, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.dpi': 600, 'figure.facecolor': 'white', 'axes.facecolor': 'white'})

def xy(row):
    return (float(row['nominal_PEA_pct']) / 100 + float(row['nominal_FA_pct']) / 200, HEIGHT * float(row['nominal_FA_pct']) / 100)

def triangle(ax):
    for fraction in (0.2, 0.4, 0.6, 0.8):
        ax.plot([fraction / 2, 1 - fraction / 2], [HEIGHT * fraction] * 2, color='#e5e5e5', lw=0.65, zorder=0)
        ax.plot([fraction, (1 + fraction) / 2], [0, HEIGHT * (1 - fraction)], color='#e5e5e5', lw=0.65, zorder=0)
        ax.plot([1 - fraction, (1 - fraction) / 2], [0, HEIGHT * (1 - fraction)], color='#e5e5e5', lw=0.65, zorder=0)
    ax.plot([0, 0.5, 1, 0], [0, HEIGHT, 0, 0], color='black', lw=1.25, zorder=1)
    ax.text(0.395, HEIGHT + 0.015, 'FAPbI$_3$', ha='right', va='center', fontsize=13, weight='bold')
    ax.text(-0.02, -0.082, 'BDAPbI$_4$', ha='center', va='top', fontsize=13, weight='bold')
    ax.text(1.025, -0.082, 'PEA$_2$PbI$_4$', ha='center', va='top', fontsize=13, weight='bold')
    ax.set(xlim=(-0.13, 1.13), ylim=(-0.11, HEIGHT + 0.1), aspect='equal')
    ax.axis('off')

def pl_triangle(final, matched, retention=False):
    fig = plt.figure(figsize=(6.0, 4.85))
    ax = fig.add_axes([0.06, 0.14, 0.75, 0.84])
    triangle(ax)
    field = 'retention_I22_over_series_max' if retention else 'acquisition22_valid_observed_max_counts'
    values = pd.to_numeric(final[field], errors='raise').to_numpy(float)
    assert np.isfinite(values).all(), 'Invalid wells require explicit categorical handling; do not silently omit them'
    norm = colors.Normalize(0, 1 if retention else float(values.max()))
    cmap = plt.get_cmap('jet')
    for repeated_column_12 in (True, False):
        rows = final.loc[final.display_marker_part.eq('right') == repeated_column_12]
        locations = np.array([xy(row) for _, row in rows.iterrows()])
        ax.scatter(*locations.T, s=24, c=cmap(norm(rows[field].astype(float))), edgecolors='none', zorder=3)
    indexed = final.set_index('well_id', drop=False)
    for well in matched.well_id:
        x, y = xy(indexed.loc[well])
        ax.scatter([x], [y], s=55, facecolors='none', edgecolors='black', linewidths=0.85, zorder=4)
    offsets = {'G1': (0.025, 0.017), 'B11': (0.02, 0.02), 'C12': (0.013, 0.023), 'H11': (0.026, -0.003), 'D5': (0.018, 0.027)}
    for number, well in enumerate(LABEL_WELLS, 1):
        x, y = xy(indexed.loc[well])
        dx, dy = offsets[well]
        ax.text(x + dx, y + dy, str(number), fontsize=10, weight='bold', zorder=5)
    color_axis = fig.add_axes([0.845, 0.19, 0.025, 0.72])
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=color_axis)
    cb.solids.set_rasterized(False)
    cb.ax.yaxis.set_major_locator(MaxNLocator(5))
    cb.set_label('PL retention at acquisition 22, $I_{22}/I_{\\max}$' if retention else 'Top-read PL peak intensity at acquisition 22 (counts)', fontsize=9, labelpad=8)
    handle = Line2D([], [], ls='none', marker='o', mfc='white', mec='black', ms=6, label='Matched photoKPFM subset')
    fig.legend(handles=[handle], loc='lower center', bbox_to_anchor=(0.42, 0.025), frameon=True, fontsize=9, borderpad=0.4)
    return fig

def normalized_traces(series, final):
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    fig.subplots_adjust(left=0.13, right=0.72, top=0.96, bottom=0.155)
    lookup = final.set_index('well_id')
    for well in LABEL_WELLS:
        part = series[series.well_id.eq(well)].sort_values('acquisition_number')
        values = pd.to_numeric(part.normalized_peak_over_all22_max, errors='raise')
        assert len(part) == 22 and np.isfinite(values).all()
        row = lookup.loc[well]
        label = '/'.join((str(int(row[k])) for k in ('nominal_FA_pct', 'nominal_BDA_pct', 'nominal_PEA_pct')))
        ax.plot(part.acquisition_number, values, color=TRACE_COLORS[well], marker='o', markersize=3.5, lw=1.5, label=label)
    ax.axvline(22, color='#555555', ls='--', lw=1, label='Acquisition 22')
    ax.set(xlim=(0.7, 22.9), ylim=(0, 1.045), xticks=[1, 5, 10, 15, 20, 22], xlabel='PL acquisition number', ylabel='Normalized PL peak intensity')
    ax.grid(color='#e5e5e5', linewidth=0.7)
    ax.set_axisbelow(True)
    ax.legend(title='Nominal composition\n(FA / BDA / PEA, %)', loc='upper left', bbox_to_anchor=(1.015, 1.0), frameon=False, fontsize=10, title_fontsize=10)
    return fig

def pl_marker_area(fraction):
    return 36 + 180 * fraction

def crossmodal_triangle(matched):
    fig = plt.figure(figsize=(6.0, 5.05))
    ax = fig.add_axes([0.07, 0.05, 0.73, 0.88])
    triangle(ax)
    values = matched.AFM_Figure5_abs_mean_Gaussian_V.to_numpy(float)
    norm = colors.Normalize(0, float(values.max()))
    cmap = plt.get_cmap('jet')
    for _, row in matched.iterrows():
        x, y = xy(row)
        ax.scatter([x], [y], s=pl_marker_area(float(row.PL_norm_to_matched13_max)), c=[cmap(norm(float(row.AFM_Figure5_abs_mean_Gaussian_V)))], edgecolors='black', lw=0.75, zorder=3)
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=fig.add_axes([0.84, 0.145, 0.025, 0.64]))
    cb.solids.set_rasterized(False)
    cb.ax.yaxis.set_major_locator(MaxNLocator(5))
    cb.set_label('Relative $|\\Delta SP|$ (V)', fontsize=10, labelpad=8)
    sizes = [0.1, 0.5, 1.0]
    handles = [Line2D([], [], ls='none', marker='o', markersize=np.sqrt(pl_marker_area(v)), markerfacecolor='#bbbbbb', markeredgecolor='black', markeredgewidth=0.6) for v in sizes]
    fig.legend(handles, [str(v) for v in sizes], title='Normalized PL peak intensity ($I_{PL}/I_{PL,\\max}$)', ncol=3, loc='upper center', bbox_to_anchor=(0.42, 0.95), frameon=False, fontsize=9, title_fontsize=10, borderpad=0.2, labelspacing=0.55, columnspacing=1.3, handletextpad=0.5)
    return fig

def rank_color(pl_rank, afm_rank):
    if pl_rank <= 4 and afm_rank <= 4:
        return '#2a9d8f'
    if afm_rank - pl_rank >= 5:
        return '#e9c46a'
    if pl_rank - afm_rank >= 5:
        return '#e76f51'
    return '#b2b2b2'

def rank_dumbbell(matched):
    rows = matched.sort_values(['PL_rank', 'well_id'], kind='stable')
    fig, ax = plt.subplots(figsize=(8.0, 5.15))
    fig.subplots_adjust(left=0.21, right=0.97, top=0.89, bottom=0.23)
    y = np.arange(len(rows))
    for j, (_, row) in enumerate(rows.iterrows()):
        if j % 2 == 0:
            ax.axhspan(j - 0.5, j + 0.5, facecolor='#f3f5f6', zorder=0)
        ax.plot([row.PL_rank, row.AFM_rank], [j, j], color=rank_color(row.PL_rank, row.AFM_rank), lw=2.2, zorder=2)
    ax.scatter(rows.PL_rank, y, marker='s', color='#4c78a8', edgecolors='white', linewidths=0.5, s=33, zorder=3, label='PL peak intensity')
    ax.scatter(rows.AFM_rank, y, marker='o', color='#f28e2b', edgecolors='white', linewidths=0.5, s=35, zorder=3, label='photoKPFM $|\\Delta SP|$')
    labels = [' / '.join((f'{int(row[k]):2d}' for k in ['nominal_FA_pct', 'nominal_BDA_pct', 'nominal_PEA_pct'])) for _, row in rows.iterrows()]
    ax.set(yticks=y, yticklabels=labels, ylim=(len(rows) - 0.4, -0.6), xlim=(0.6, 13.5), xticks=[1, 3, 5, 7, 9, 11, 13], xlabel='Response rank among 13 nominal matches (1 = strongest)')
    ax.axvline(4.5, color='#77bfc3', lw=0.7, ls='--')
    ax.tick_params(axis='y', length=0, labelsize=9)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.grid(axis='x', color='#e2e5e8', lw=0.7)
    ax.set_axisbelow(True)
    ax.text(-0.23, 1.035, 'FA / BDA / PEA (%)', transform=ax.transAxes, ha='left', fontsize=9, weight='bold')
    ax.legend(loc='lower left', bbox_to_anchor=(-0.015, 1.02), ncol=2, frameon=False, fontsize=9)
    handles = [Line2D([], [], color=color, lw=2.2, label=label) for color, label in [('#2a9d8f', 'Both in top 4'), ('#e9c46a', 'PL leads by ≥ 5 ranks'), ('#e76f51', 'photoKPFM leads by ≥ 5 ranks'), ('#b2b2b2', 'Other')]]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.59, 0.045), ncol=2, fontsize=8, frameon=False)
    return fig

def crossmodal_scatter(matched):
    fig, ax = plt.subplots(figsize=(7.1, 4.3))
    fig.subplots_adjust(left=0.13, right=0.97, top=0.95, bottom=0.16)
    ax.errorbar(matched.PL_raw_max_counts / 1000, matched.AFM_S5_mean_abs_Gaussian_V, yerr=matched.AFM_archived_delta_SP_se_V, fmt='o', ms=np.sqrt(43), mfc='#4c78a8', mec='black', mew=0.65, ecolor='#999999', elinewidth=0.85, capsize=3, zorder=3)
    ax.set(xlabel='Top-read PL peak intensity at acquisition 22 (×10³ counts)', ylabel='Relative $|\\Delta SP|$ (V)')
    ax.grid(color='#e5e5e5', lw=0.7)
    ax.set_axisbelow(True)
    ax.margins(0.06)
    return fig

def render_all(final, series, matched, output):
    style()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    makers = {"4a": lambda: pl_triangle(final, matched), "4b": lambda: pl_triangle(final, matched, True), "S4": lambda: normalized_traces(series, final), "5a": lambda: crossmodal_triangle(matched), "5b": lambda: rank_dumbbell(matched), "S5": lambda: crossmodal_scatter(matched)}
    outputs = []
    for panel, make in makers.items():
        figure = make()
        for extension in ("svg", "pdf", "png"):
            destination = output / ("Figure" + panel + "_gain60." + extension)
            figure.savefig(destination, dpi=600, facecolor="white", transparent=False)
            outputs.append(destination)
        plt.close(figure)
    return outputs
