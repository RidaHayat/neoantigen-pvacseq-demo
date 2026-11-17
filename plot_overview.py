import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

df = pd.read_csv("/pvacseq_output_data/top_candidates.csv")
mt = df["Best MT IC50 Score"]
wt = df["Corresponding WT IC50 Score"]
fc = df["Corresponding Fold Change"]
expr = df["Gene Expression"]

fig, (a, b) = plt.subplots(1, 2, figsize=(14, 5.5))

# Panel A: binding vs fold change, bubble size = expression
a.scatter(mt, fc, s=30 + expr.clip(upper=400) / 2, alpha=0.6)
for g, x, y in zip(df["Gene Name"], mt, fc):
    a.annotate(g, (x, y), textcoords="offset points", xytext=(5, 5), fontsize=8)
a.set_xscale("log")
a.set_yscale("log")
a.axvline(500, ls="--", c="grey")
a.axhline(1, ls=":", c="grey")
a.set_xlabel("Predicted mutant IC50 (nM), lower = stronger binding")
a.set_ylabel("Fold change (normal IC50 / mutant IC50)")
a.set_title("A. Binding strength vs. difference from normal")

# Panel B: mutant vs normal peptide binding per gene
d = df.assign(wt=wt, mt=mt, fc=fc).sort_values("fc")
for i, (g, m, w) in enumerate(zip(d["Gene Name"], d["mt"], d["wt"])):
    b.plot([m, w], [i, i], c="lightgrey", zorder=1)
b.scatter(d["mt"], range(len(d)), label="Mutant peptide", zorder=2)
b.scatter(d["wt"], range(len(d)), marker="s", label="Normal peptide", zorder=2)
b.set_yticks(range(len(d)))
b.set_yticklabels(d["Gene Name"])
b.set_xscale("log")
b.axvline(500, ls="--", c="grey")
b.set_xlabel("Predicted IC50 (nM)")
b.set_title("B. Mutant vs. normal peptide binding")
b.legend()

plt.tight_layout()
plt.savefig("/pvacseq_output_data/neoantigen_overview.png", dpi=150)
print("Saved neoantigen_overview.png")
