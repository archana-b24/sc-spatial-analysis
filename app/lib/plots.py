"""Reusable figures. UMAPs use WebGL (scattergl) on the stratified subsample stored by `export`."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .ui import CONDITION_COLOR, CONDITION_LABEL, GREY, SERIES, style

COMPARTMENT_COLOR = {
    "Epithelial": SERIES[0],
    "Immune": SERIES[1],
    "Endothelial": SERIES[2],
    "Stroma": SERIES[3],
    "Unknown": GREY,
}


def umap_categorical(
    df: pd.DataFrame, col: str, colors: dict | None = None, highlight: str | None = None, height: int = 560
) -> go.Figure:
    """Categorical UMAP. With > 8 categories (cell types, studies) colour cannot carry identity safely, so
    the caller passes `highlight` and everything else is grey."""
    fig = go.Figure()
    if highlight is not None:
        groups = [
            ("other", df[df[col] != highlight], GREY, 3),
            (highlight, df[df[col] == highlight], SERIES[0], 4),
        ]
    else:
        cats = [c for c in (colors or {}) if c in set(df[col])] or sorted(df[col].unique())
        groups = [(c, df[df[col] == c], (colors or {}).get(c, GREY), 3) for c in cats]
    for name, d, color, size in groups:
        label = CONDITION_LABEL.get(name, name)
        fig.add_trace(
            go.Scattergl(
                x=d["umap_1"],
                y=d["umap_2"],
                mode="markers",
                name=f"{label} ({len(d):,})",
                marker=dict(color=color, size=size, opacity=0.75),
                customdata=np.stack([d["harmonised_label"], d["condition"], d["study"]], axis=1),
                hovertemplate="%{customdata[0]}<br>%{customdata[1]} · %{customdata[2]}<extra></extra>",
            )
        )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False, scaleanchor="x")
    return style(fig, height)


def umap_continuous(df: pd.DataFrame, values: pd.Series, title: str, height: int = 560) -> go.Figure:
    order = np.argsort(values.to_numpy())  # high values drawn last (on top)
    d, v = df.iloc[order], values.iloc[order]
    fig = go.Figure(
        go.Scattergl(
            x=d["umap_1"],
            y=d["umap_2"],
            mode="markers",
            marker=dict(
                color=v,
                size=3,
                colorscale=[[0, "#e8e7e3"], [0.25, "#86b6ef"], [1, "#0d366b"]],
                colorbar=dict(title=title, thickness=12),
            ),
            customdata=np.stack([d["harmonised_label"], d["condition"], v.round(2)], axis=1),
            hovertemplate="%{customdata[0]} · %{customdata[1]}<br>"
            + title
            + " %{customdata[2]}<extra></extra>",
            showlegend=False,
        )
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False, scaleanchor="x")
    return style(fig, height)


def condition_bar(df: pd.DataFrame, x: str, y: str, orientation: str = "h", height: int = 360) -> go.Figure:
    """Grouped bars, one colour per condition (same colours on every page)."""
    fig = go.Figure()
    for cond, color in CONDITION_COLOR.items():
        d = df[df["condition"] == cond]
        if d.empty:
            continue
        kw = dict(x=d[x], y=d[y]) if orientation == "v" else dict(x=d[y], y=d[x])
        fig.add_trace(go.Bar(name=CONDITION_LABEL[cond], marker_color=color, orientation=orientation, **kw))
    fig.update_layout(barmode="group", bargap=0.25, bargroupgap=0.08)
    return style(fig, height)


def heatmap(mat: pd.DataFrame, title: str = "", diverging: bool = True, height: int = 480, fmt: str = ".2f"):
    from .ui import DIVERGING

    z = mat.to_numpy(dtype=float)
    kw = {}
    if diverging:
        m = np.nanmax(np.abs(z)) or 1
        kw = dict(colorscale=DIVERGING, zmin=-m, zmax=m)
    else:
        kw = dict(colorscale=[[0, "#f0efec"], [1, "#1c5cab"]])
    fig = go.Figure(
        go.Heatmap(
            z=z,
            x=list(mat.columns),
            y=list(mat.index),
            colorbar=dict(title=title, thickness=12),
            hovertemplate="%{y} · %{x}<br>" + title + " %{z:" + fmt + "}<extra></extra>",
            xgap=1,
            ygap=1,
            **kw,
        )
    )
    fig.update_xaxes(tickangle=-45)
    return style(fig, height)
