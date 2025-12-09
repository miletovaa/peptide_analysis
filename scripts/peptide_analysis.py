import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import probplot, shapiro, f_oneway, kruskal, gaussian_kde
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA

from utils.config import cfg
from utils.output import create_output_dir, get_output_file_name, save_figure

pd.set_option('future.no_silent_downcasting', True)


OUTPUT_DIR = create_output_dir(cfg)


### Loading file to dataframe
with pd.ExcelFile(os.path.join(cfg.get("FILE_PATH_INPUT"), cfg.data.file_name_input)) as xls:
    sheet_names = xls.sheet_names
    final_sheets = [name for name in sheet_names if name.upper().startswith("FINAL")]

    if not final_sheets:
        raise ValueError(f"No sheet starting with 'FINAL' found in {cfg.data.file_name_input}. Available sheets: {sheet_names}")
    
    df = pd.read_excel(xls, sheet_name=final_sheets[0])
    print(f"Loaded sheet '{final_sheets[0]}' from {cfg.data.file_name_input}")

# Fixing double-row header
df = df[1:]

# %% [markdown]
# ### Step 1: Dataframe cleaning and preview

# %%
# Forward fill missing values (aka fixing merged cells with protein names)
df["Protein"] = df["Protein"].ffill()
print(f"Successfully performed forward fill on missing values.")

print("Total number of peptides before filtering:", len(df))

sample_prefix_arr = ("BIOPD", "RBD", "CON")
value_cols = [c for c in df.columns if c.startswith(sample_prefix_arr)]

df[value_cols] = df[value_cols].replace("<LOQ", 0.0).astype(float)

groups = {
    "PD": [c for c in value_cols if c.startswith("BIOPD")],
    "RBD": [c for c in value_cols if c.startswith("RBD")],
    "CON": [c for c in value_cols if c.startswith("CON")]
}

print(f"Step 1: Preview of the dataframe:")
df


# %% [markdown]
# ### Step 2: Removing samples with more than specified threshold of zero values (<LOQ; default: 60%)

# %%
zero_threshold = cfg.preprocessing.zero_threshold if cfg.preprocessing.zero_threshold is not None else 0.6

zero_fractions = (df.eq(0.0) | df.isna()).sum(axis=1) / len(df)

df_filtered = df[zero_fractions <= zero_threshold]
print(f"Step 2: Filtered peptides that have more than %d%% of '<LOQ': kept %d of %d" % (zero_threshold * 100, len(df_filtered), len(df)))
df = df_filtered.copy()


# %% [markdown]
# ### Step 3: Imputing remaining '<LOQ' values with LOQ/2 for each peptide

# %%
for col in value_cols:
    min_detected = df[col].min(skipna=True)
    impute_value = min_detected / 2
    df.fillna({col: impute_value}, inplace=True)

print(f"Step 3: Imputed remaining '<LOQ' values with LOQ/2 for each peptide.")

# save intermediate result
if cfg.run.save_intermediate:
    intermediate_file = get_output_file_name(step="s03", type="prep", description="zero_thresholded_imputed", ext="csv", cfg=cfg, output_dir=OUTPUT_DIR)
    df.to_csv(intermediate_file, index=False)
    print(f"Intermediate result, zero-thresholded and imputed, saved to {intermediate_file}")

# %% [markdown]
# ### Step 4: Shapiro-Wilk test for normality on each peptide across samples

# %%
def do_shapiro(df, groups):
    results = []
    for _, row in df.iterrows():
        for group, cols in groups.items():
            values = row[cols].astype(float).values
            if len(values) < 3:
                continue

            stat, pval = shapiro(values)
            _, (slope, intercept, r) = probplot(values)
            
            results.append({
                "peptide": row["Peptide"],
                "protein": row["Protein"],
                "group": group,

                "is_normal_shapiro": pval > 0.01,
                "is_normal_qq": r >= 0.990,

                "shapiro_stat": stat,
                "shapiro_pval": pval,

                "qq_slope": slope,
                "qq_intercept": intercept,
                "qq_rvalue": r
            })
    return results

# * Applied Shapiro-Wilk test on INITIAL peptide data (before log transformation)
# ? normality_df: peptide, protein, group, is_normal_shapiro, is_normal_qq, 
# ?               shapiro_stat, shapiro_pval, qq_slope, qq_intercept, qq_rvalue
normality_df = pd.DataFrame(do_shapiro(df, groups))
normality_df.sort_values(by=["peptide", "group"], inplace=True) # by=["group", "peptide"]

# * Count percentage of non-normal peptides
# ? normal: array of peptides and their groups that have is_normal_shapiro == True
normal = normality_df.loc[normality_df["is_normal_shapiro"], ["peptide", "group"]]
print(f"Number of normal group-peptides: %d out of %d (%.2f%%)" % (normal.shape[0], normality_df.shape[0], normal.shape[0] / normality_df.shape[0] * 100))

# * Peptides normal in all 3 groups (by Shapiro)
normal_all_3 = (
    normality_df.groupby("peptide")
    .agg(n_groups=("group", "nunique"), all_normal=("is_normal_shapiro", "all"))
    .query("n_groups == 3 and all_normal")
    .index.tolist()
)
print(f"Peptides normal in all 3 groups (Shapiro): %d" % len(normal_all_3))


# * Apply log10 transformation and re-test normality
df_log = df.copy()
df_log[value_cols] = np.log10(
    df[value_cols].astype(float).clip(lower=1e-12)
)
normality_df_log = pd.DataFrame(do_shapiro(df_log, groups))
normality_df_log.sort_values(by=["peptide", "group"], inplace=True)

# * Count percentage of non-normal peptides after log transformation
normal_log = normality_df_log.loc[normality_df_log["is_normal_shapiro"], ["peptide", "group"]]
print(f"Number of normal group-peptides (after log transformation): %d out of %d (%.2f%%)" % (normal_log.shape[0], normality_df_log.shape[0], normal_log.shape[0] / normality_df_log.shape[0] * 100))

# * Peptides normal in all 3 groups (by Shapiro)
normal_all_3 = (
    normality_df_log.groupby("peptide")
    .agg(n_groups=("group", "nunique"), all_normal=("is_normal_shapiro", "all"))
    .query("n_groups == 3 and all_normal")
    .index.tolist()
)
print(f"Peptides normal in all 3 groups (Shapiro): %d" % len(normal_all_3))

print(f"Step 4: Completed Shapiro-Wilk tests for each peptide within each group.")

if cfg.run.save_intermediate:
    normality_df.to_csv(
        get_output_file_name(step="s04", type="stats", description="shapiro_normality", ext="csv", cfg=cfg, output_dir=OUTPUT_DIR),
        index=False
    )
    normality_df_log.to_csv(
        get_output_file_name(step="s04", type="stats", description="shapiro_normality_log_transformed", ext="csv", cfg=cfg, output_dir=OUTPUT_DIR),
        index=False
    )
    print("Saved normality test results to CSV files.")

# %% [markdown]
# ### Step 5: Doing statistical tests between groups

# %%
# TODO: condition for choosing test based on normality results
# * Using Kruskal-Wallis test for all peptides (non-parametric)
kw_results = []
group_names = list(groups.keys())
for _, row in df.iterrows():
    values_by_group = [row[groups[g]].astype(float).values for g in group_names]
    stat, pval = kruskal(*values_by_group)
    kw_results.append({
        "peptide": row["Peptide"],
        "protein": row["Protein"],
        "kw_stat": stat,
        "kw_pval": pval,
        "is_significant": pval < 0.01,
    })

kw_df = pd.DataFrame(kw_results)
kw_df.sort_values(by=["is_significant", "kw_pval", "peptide"], ascending=[False, True, True], inplace=True)


# * Using ANOVA for peptides normal in all 3 groups
anova_results = []
for peptide in normal_all_3:
    row = df.loc[df["Peptide"] == peptide].iloc[0]
    values_by_group = [row[groups[g]].astype(float).values for g in group_names]
    stat, pval = f_oneway(*values_by_group)
    anova_results.append({
        "peptide": row["Peptide"],
        "protein": row["Protein"],
        "anova_stat": stat,
        "anova_pval": pval,
        "is_significant": pval < 0.01,
    })
anova_df = pd.DataFrame(anova_results)
anova_df.sort_values(by=["is_significant", "anova_pval", "peptide"], ascending=[False, True, True], inplace=True)

print(f"Step 5: Statistical tests between groups done.")

if cfg.run.save_intermediate:
    kw_df.to_csv(
        get_output_file_name(step="s05", type="stats", description="kruskal_wallis_results", ext="csv", cfg=cfg, output_dir=OUTPUT_DIR),
        index=False
    )
    print("Saved statistical test results to CSV files.")


# %% [markdown]
# ### Step 6: Applying Bonferroni and Benjamini-Hochberg correction and detecting significant peptides
# 

# %%
sig_kw = set(kw_df.loc[kw_df["is_significant"], "peptide"])
sig_anova = set(anova_df.loc[anova_df["is_significant"], "peptide"])

sig_peptides = sorted(sig_kw.union(sig_anova))
n_sig = len(sig_peptides)

kw_pvals = kw_df["kw_pval"].values

_, kw_pvals_bonf, _, _ = multipletests(kw_pvals, method="bonferroni")
_, kw_pvals_bh, _, _    = multipletests(kw_pvals, method="fdr_bh")

kw_df["kw_pval_bonf"] = kw_pvals_bonf
kw_df["kw_pval_fdr_bh"] = kw_pvals_bh

kw_df["is_significant_bonf"] = kw_df["kw_pval_bonf"] < 0.01
kw_df["is_significant_fdr"]  = kw_df["kw_pval_fdr_bh"] < 0.01
sig_kw_fdr = set(kw_df.loc[kw_df["is_significant_fdr"], "peptide"])
sig_peptides_fdr = sorted(sig_kw_fdr)

sig_kw_bonf = set(kw_df.loc[kw_df["is_significant_bonf"], "peptide"])
sig_peptides_bonf = sorted(sig_kw_bonf)
if n_sig == 0:
    print(f"Step 6: No significant peptides found. Skipping visualizations.")

print(f"Step 6: Corrections on the statistics done.")
print(f"Unique significant peptides from Kruskal/ANOVA after Bonferroni correction:", len(sig_peptides_bonf))
print(f"Unique significant peptides from Kruskal/ANOVA after FDR correction:", len(sig_peptides_fdr))

# if cfg.run.use_correction:
#     sig_peptides = sig_peptides_bonf

print(sig_peptides)
print(len(sig_peptides))

if cfg.run.save_intermediate:
    kw_df.to_csv(
        get_output_file_name(step="s06", type="stats", description="kruskal_wallis_results_with_corrections", ext="csv", cfg=cfg, output_dir=OUTPUT_DIR),
        index=False
    )
    print("Saved Kruskal-Wallis results with corrections to CSV file.")

# %% [markdown]
# ### Step 7: Visualizations of the groups of peptides that are significant
# 
# **7.1:** Set colors for groups, got labels for peptides with their proteins
# 

# %%
group_colors = {
    "PD":  "#636EFA",   # blue
    "RBD": "#EF553B",   # red
    "CON": "#00CC96",   # green
}

group_names = list(groups.keys())

peptide_to_protein = (
    df.loc[:, ["Peptide", "Protein"]]
    .dropna(subset=["Peptide"])
    .drop_duplicates(subset=["Peptide"], keep="first")
    .set_index("Peptide")["Protein"]
    .to_dict()
)

def _title_with_protein(peptide_name: str) -> str:
    protein = peptide_to_protein.get(peptide_name, "Unknown protein")
    return f"{protein} ({peptide_name})"

peptide_titles = [_title_with_protein(p) for p in sig_peptides]


# %% [markdown]
# **7.2**: Box plots

# %%
fig_box = make_subplots(
    rows=1,
    cols=len(sig_peptides),
    shared_yaxes=False,
    subplot_titles=peptide_titles
)
for col_idx, peptide in enumerate(sig_peptides, start=1):
    row_pep = df.loc[df["Peptide"] == peptide].iloc[0]
    for g in group_names:
        values = row_pep[groups[g]].astype(float).values
        fig_box.add_trace(
            go.Box(
                y=values,
                name=g,
                boxpoints="outliers",
                showlegend=(col_idx == 1),
                marker_color=group_colors.get(g),
            ),
            row=1,
            col=col_idx
        )
width_box = max(350 * len(sig_peptides), 800)
height_box = 400
fig_box.update_layout(
    title="Significant Peptides – Box Plots per Peptide (3 groups each)",
    width=width_box,
    height=height_box
)
# Axis labels for boxplots
fig_box.update_xaxes(title_text="Group")
fig_box.update_yaxes(title_text="Peptide intensity", row=1, col=1)
save_figure(
    fig_box,
    "s07",
    "significant_peptides_boxplots",
    cfg,
    OUTPUT_DIR
)
print(f"Step 7.2: Box plots done.")

# %% [markdown]
# **7.3**: Histograms: median bars

# %%
fig_hist = make_subplots(
    rows=1,
    cols=len(sig_peptides),
    shared_yaxes=False,
    subplot_titles=peptide_titles
)
for col_idx, peptide in enumerate(sig_peptides, start=1):
    row_pep = df.loc[df["Peptide"] == peptide].iloc[0]
    for g in group_names:
        values = row_pep[groups[g]].astype(float).values
        values = values[~np.isnan(values)]
        if len(values) == 0:
            median_val = np.nan
        else:
            median_val = float(np.median(values))
        fig_hist.add_trace(
            go.Bar(
                x=[g],
                y=[median_val],
                name=g,
                marker_color=group_colors.get(g),
                showlegend=(col_idx == 1)
            ),
            row=1,
            col=col_idx
        )
fig_hist.update_layout(
    barmode="group"
)
width_hist = max(300 * len(sig_peptides), 300)
height_hist = 350
fig_hist.update_layout(
    title="Significant Peptides – Group Medians per Peptide (3 bars each)",
    width=width_hist,
    height=height_hist
)
# Axis labels for median-bar histograms
fig_hist.update_xaxes(title_text="Group")
fig_hist.update_yaxes(title_text="Median peptide intensity", row=1, col=1)
save_figure(
    fig_hist,
    "s07",
    "significant_peptides_median_bars",
    cfg,
    OUTPUT_DIR
)
print(f"Step 7.3: Median bar plots done.")


# %% [markdown]
# **7.4**: Density distributions

# %%
fig_dist = make_subplots(
    rows=1,
    cols=len(sig_peptides),
    shared_yaxes=False,
    subplot_titles=peptide_titles
)
for col_idx, peptide in enumerate( sig_peptides, start=1):
    row_pep = df.loc[df["Peptide"] == peptide].iloc[0]

    all_vals = []
    group_vals = {}
    for g in group_names:
        vals = row_pep[groups[g]].astype(float).values
        vals = vals[~np.isnan(vals)]
        if len(vals) == 0:
            continue
        group_vals[g] = vals
        all_vals.append(vals)
    if not all_vals:
        continue
    all_vals = np.concatenate(all_vals)
    x_min, x_max = float(all_vals.min()), float(all_vals.max())
    if x_min == x_max:
        x_min -= 1e-6
        x_max += 1e-6
    xs = np.linspace(x_min, x_max, 200)

    for g in group_names:
        if g not in group_vals:
            continue
        vals = group_vals[g]
        if len(vals) < 2:
            continue
        kde = gaussian_kde(vals)
        ys = kde(xs)
        fig_dist.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                name=g,
                showlegend=(col_idx == 1),
                line=dict(color=group_colors.get(g)),
            ),
            row=1,
            col=col_idx
        )
width_dist = max(400 * len(sig_peptides), 400)
height_dist = 350
fig_dist.update_layout(
    title="Significant Peptides – Density Distributions per Peptide (3 groups each)",
    width=width_dist,
    height=height_dist
)
# Axis labels for density plots
fig_dist.update_xaxes(title_text="Peptide intensity")
fig_dist.update_yaxes(title_text="Density", row=1, col=1)
save_figure(
    fig_dist,
    "s07",
    "significant_peptides_density_distributions",
    cfg,
    OUTPUT_DIR
)
print(f"Step 7.4: Distribution plots saved.")


# %% [markdown]
# ### Step 8: Analysis of RBD subgroups of significant peptides

# %%
rbd_groups = {
    "A": [c for c in value_cols if c.endswith("RBDA")],
    "B": [c for c in value_cols if c.endswith("RBDB")],
    "C": [c for c in value_cols if c.endswith("RBDC")],
}

rbd_kw_results = []
rbd_descriptives = []

df_by_peptide = df.set_index("Peptide", drop=False)

for peptide in df_by_peptide.index.unique():
    row_pep = df_by_peptide.loc[peptide]
    if isinstance(row_pep, pd.DataFrame):
        row_pep = row_pep.iloc[0]

    subgroup_values = []
    for subgroup, cols in rbd_groups.items():
        vals = row_pep[cols].astype(float).values
        vals = vals[~np.isnan(vals)]
        subgroup_values.append(vals)
        rbd_descriptives.append({
            "peptide": peptide,
            "rbd_group": subgroup,
            "n": len(vals),
            "median": float(np.median(vals)),
            "mean": float(np.mean(vals)),
        })

    stat, pval = kruskal(*subgroup_values)
    rbd_kw_results.append({
        "peptide": peptide,
        "kw_stat": stat,
        "kw_pval": pval,
        "is_significant": pval < 0.01,
    })

    rbd_kw_df = pd.DataFrame(rbd_kw_results)
    rbd_kw_df.sort_values(by=["is_significant", "kw_pval", "peptide"], ascending=[False, True, True], inplace=True)

print("Step 8: RBD subgroup analysis completed. %d peptides tested." % (len(rbd_kw_results)))

print(rbd_kw_df)
if cfg.run.save_intermediate:
    rbd_kw_df.to_csv(
        get_output_file_name(step="s08", type="stats", description="rbd_subgroup_kruskal_wallis_results", ext="csv", cfg=cfg, output_dir=OUTPUT_DIR),
        index=False
    )
    rbd_desc_df = pd.DataFrame(rbd_descriptives)
    rbd_desc_df.to_csv(
        get_output_file_name(step="s08", type="stats", description="rbd_subgroup_descriptives", ext="csv", cfg=cfg, output_dir=OUTPUT_DIR),
        index=False
    )
    print("Saved RBD subgroup analysis results to CSV files.")

# %%
print(rbd_kw_df[rbd_kw_df["is_significant"]])

# %%
print(rbd_kw_df.loc[rbd_kw_df["peptide"].isin(sig_peptides)])

# %% [markdown]
# # PCA

# %%
try:
    data_matrix = df[value_cols].astype(float).T
    sample_ids = data_matrix.index.to_list()
    data_log = np.log10(data_matrix + 1.0e-6)
    if "Peptide" in df.columns:
        raw_labels = df["Peptide"].astype(str).fillna("Unknown").tolist()
    else:
        raw_labels = [f"feature_{i}" for i in range(data_matrix.shape[1])]
    if len(raw_labels) != data_matrix.shape[1]:
        raw_labels = [f"feature_{i}" for i in range(data_matrix.shape[1])]
    data_matrix.columns = raw_labels

    # Pareto scaling
    col_means = data_log.mean(axis=0)
    col_std = data_log.std(axis=0, ddof=1).replace(0, np.nan)
    pareto_scale = np.sqrt(col_std).replace(0, np.nan)
    data_pareto = (data_log - col_means) / pareto_scale
    data_pareto = data_pareto.dropna(axis=1, how="any")
    retained_features = data_pareto.columns.to_list()

    # Number of components is set here
    pca = PCA(n_components=0.95) # automatically calculate all components
    pca_scores = pca.fit_transform(data_pareto.values)

    sample_groups = []
    for sid in sample_ids:
        if sid.startswith("BIOPD"): sample_groups.append("PD")
        elif sid.startswith("RBD"): sample_groups.append("RBD")
        elif sid.startswith("CON"): sample_groups.append("CON")
        else: sample_groups.append("Other")

    pca_df = pd.DataFrame({
        "sample": sample_ids,
        "group": sample_groups,
        **{f"PC{i+1}": pca_scores[:, i] for i in range(pca.n_components_)}
    })
    print(pca_df.head())

    component_labels = [f"PC{i+1}" for i in range(pca.n_components_)]
    loadings_df = pd.DataFrame(
        pca.components_.T,
        columns=component_labels,
        index=retained_features,
    )

    # loadings_export = loadings_df.reset_index().rename(columns={"index": "feature"})
    # loadings_export.to_csv(os.path.join(OUTPUT_DIR, "pca_loadings.csv"), index=False)
    loadings_df["loading_magnitude"] = np.sqrt(loadings_df["PC1"] ** 2 + loadings_df["PC2"] ** 2)
    top_n_loadings = min(20, len(loadings_df))
    top_loadings = loadings_df.nlargest(top_n_loadings, "loading_magnitude") if top_n_loadings else loadings_df

    print("PCA performed successfully.")

    if cfg.run.save_intermediate:
        pca_df.to_csv(
            get_output_file_name(step="s09", type="viz", description="pca_scores", ext="csv", cfg=cfg, output_dir=OUTPUT_DIR),
            index=False
        )
        loadings_df.to_csv(
            get_output_file_name(step="s09", type="viz", description="pca_loadings", ext="csv", cfg=cfg, output_dir=OUTPUT_DIR),
            index=True
        )
        print("Saved PCA scores and loadings to CSV files.")
except Exception as e:
    print(f"Error during PCA computation: {e}")


# %% [markdown]
# ### PCA scores plot

# %%
try:
    group_colors_pca = {
        "PD":  "#636EFA",
        "RBD": "#EF553B",
        "CON": "#00CC96",
        "Other": "#AB63FA",
    }

    fig_pca = go.Figure()
    for grp in sorted(pca_df["group"].unique()):
        grp_df = pca_df[pca_df["group"] == grp]
        fig_pca.add_trace(
            go.Scatter(
                x=grp_df["PC1"],
                y=grp_df["PC2"],
                mode="markers",
                name=grp,
                marker=dict(
                    size=10,
                    opacity=0.8,
                    color=group_colors_pca.get(grp, "#888888"),
                ),
                text=grp_df["sample"],
                hovertemplate="Sample: %{text}<br>PC1: %{x:.3f}<br>PC2: %{y:.3f}<extra></extra>",
            )
        )

    fig_pca.update_layout(
        title="PCA Score scatter",
        xaxis_title=f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)",
        yaxis_title=f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)",
        width=800,
        height=700,
    )
    save_figure(
        fig_pca,
        "s09",
        "pca_scores_scatter",
        cfg,
        OUTPUT_DIR
    )
    print(f"Step 9: PCA plotting done.")
except Exception as e:
    print(f"Error during PCA computation or plotting: {e}")

# %% [markdown]
# ## PCA loadings plot

# %%
try:
    if not loadings_df.empty:
        fig_loadings = go.Figure()
        for feature, row in top_loadings.iterrows():
            fig_loadings.add_shape(
                type="line",
                x0=0,
                y0=0,
                x1=row["PC1"],
                y1=row["PC2"],
                line=dict(color="#BBBBBB", width=1),
            )
        fig_loadings.add_trace(
            go.Scatter(
                x=top_loadings["PC1"],
                y=top_loadings["PC2"],
                mode="markers+text",
                text=top_loadings.index,
                textposition="top center",
                marker=dict(size=9, color="#333333"),
                hovertemplate="Feature: %{text}<br>PC1 loading: %{x:.3f}<br>PC2 loading: %{y:.3f}<extra></extra>",
                showlegend=False,
            )
        )
        fig_loadings.update_layout(
            title="PCA Loadings (PC1 vs PC2)",
            xaxis_title="PC1 loading",
            yaxis_title="PC2 loading",
            xaxis=dict(zeroline=True, zerolinewidth=1, zerolinecolor="#333"),
            yaxis=dict(zeroline=True, zerolinewidth=1, zerolinecolor="#333"),
            width=800,
            height=700,
        )
        save_figure(
            fig_loadings,
            "s09",
            "pca_loadings",
            cfg,
            OUTPUT_DIR
        )
        print(f"Step 9: PCA loadings plotting done.")
except Exception as e:
    print(f"Error during PCA loadings plotting: {e}")

# %% [markdown]
# ## PCA biplot

# %%
try:
    score_max = float(np.max(np.abs(pca_scores[:, :2]))) if pca_scores.size else 1.0
    loading_max = float(np.max(np.abs(top_loadings[["PC1", "PC2"]].values))) if not top_loadings.empty else 1.0
    loading_max = loading_max if loading_max != 0 else 1.0
    biplot_scale = (score_max / loading_max) * 0.7 if loading_max else 1.0
    scaled_loadings = top_loadings[["PC1", "PC2"]] * biplot_scale

    fig_biplot = go.Figure()
    for grp in sorted(pca_df["group"].unique()):
        grp_df = pca_df[pca_df["group"] == grp]
        fig_biplot.add_trace(
            go.Scatter(
                x=grp_df["PC1"],
                y=grp_df["PC2"],
                mode="markers",
                name=f"{grp} samples",
                marker=dict(
                    size=10,
                    opacity=0.8,
                    color=group_colors_pca.get(grp, "#888888"),
                ),
                text=grp_df["sample"],
                hovertemplate="Sample: %{text}<br>PC1: %{x:.3f}<br>PC2: %{y:.3f}<extra></extra>",
            )
        )

    for feature, row in scaled_loadings.iterrows():
        fig_biplot.add_shape(
            type="line",
            x0=0,
            y0=0,
            x1=row["PC1"],
            y1=row["PC2"],
            line=dict(color="#444444", width=1),
        )

    fig_biplot.add_trace(
        go.Scatter(
            x=scaled_loadings["PC1"],
            y=scaled_loadings["PC2"],
            mode="markers+text",
            text=scaled_loadings.index,
            textposition="top center",
            marker=dict(color="#444444", size=9, symbol="diamond"),
            name="Feature loadings",
            hovertemplate="Feature: %{text}<br>Scaled PC1 loading: %{x:.3f}<br>Scaled PC2 loading: %{y:.3f}<extra></extra>",
        )
    )

    fig_biplot.update_layout(
        title="PCA Biplot (PC1 vs PC2)",
        xaxis_title=f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)",
        yaxis_title=f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)",
        width=900,
        height=750,
    )
    save_figure(
        fig_biplot,
        "s09",
        "pca_biplot",
        cfg,
        OUTPUT_DIR
    )
except Exception as e:
    print(f"Error during PCA biplot plotting: {e}")

print("Step 9: PCA plots and scores saved successfully.")

# %%
try:
    use_only_significant = False

    if use_only_significant:
        peptide_subset = sig_peptides
    else:
        peptide_subset = df["Peptide"].dropna().unique().tolist()

    pca_df_raw = df[df["Peptide"].isin(peptide_subset)]
    data_matrix = pca_df_raw.set_index("Peptide")[value_cols].astype(float)

    data_log = np.log10(data_matrix + 1e-6)
    data_log = data_log.replace([np.inf, -np.inf], np.nan)

    data_log = data_log.dropna(axis=0, how="all")

    row_means = data_log.mean(axis=1)
    data_log = data_log.T.fillna(row_means).T

    # Pareto scaling
    row_std = data_log.std(axis=1, ddof=1)
    pareto_scale = np.sqrt(row_std)
    pareto_scale[pareto_scale == 0] = 1.0
    data_centered = data_log.sub(row_means, axis=0)
    data_pareto = data_centered.div(pareto_scale, axis=0).dropna(axis=0, how="any")

    # PCA (rows = peptides, cols = samples)
    # n_components = min(data_pareto.shape[0], data_pareto.shape[1])
    pca = PCA(n_components=2)
    pca_scores = pca.fit_transform(data_pareto.values)

    peptides_final = data_pareto.index.tolist()
    proteins_final = (
        df.loc[df["Peptide"].isin(peptides_final), ["Peptide", "Protein"]]
        .drop_duplicates("Peptide")
        .set_index("Peptide")["Protein"]
        .to_dict()
    )

    peptide_pca_df = pd.DataFrame({
        "peptide": peptides_final,
        "protein": [proteins_final.get(p, "Unknown") for p in peptides_final],
        "PC1": pca_scores[:, 0],
        "PC2": pca_scores[:, 1],
    })

    # Coloring by protein
    protein_colors = {}
    unique_proteins = peptide_pca_df["protein"].unique()
    for i, pr in enumerate(unique_proteins):
        protein_colors[pr] = f"hsl({(i * 47) % 360},70%,50%)"

    fig_pep_pca = go.Figure()
    for pr in unique_proteins:
        temp = peptide_pca_df[peptide_pca_df["protein"] == pr]
        fig_pep_pca.add_trace(
            go.Scatter(
                x=temp["PC1"],
                y=temp["PC2"],
                mode="markers",
                name=pr,
                marker=dict(size=8, opacity=0.8, color=protein_colors[pr]),
                text=temp["peptide"],
                hovertemplate="Peptide: %{text}<br>PC1: %{x:.3f}<br>PC2: %{y:.3f}<extra></extra>",
            )
        )

    fig_pep_pca.update_layout(
        title="Peptide-level PCA (log10 + Pareto scaled)",
        xaxis_title=f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)",
        yaxis_title=f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)",
        width=900,
        height=800,
    )
    save_figure(
        fig_pep_pca,
        "s09",
        "pca_peptide_level",
        cfg,
        OUTPUT_DIR
    )
except Exception as e:
    print(f"Error during peptide-level PCA: {e}")
