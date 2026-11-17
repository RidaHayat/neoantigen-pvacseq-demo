import pandas as pd

d = "/pvacseq_output_data/MHC_Class_I/HCC1395_TUMOR_DNA.MHC_I."
allep = pd.read_csv(d + "all_epitopes.tsv", sep="\t", low_memory=False)
agg = pd.read_csv(d + "all_epitopes.aggregated.tsv", sep="\t")
filt = pd.read_csv(d + "filtered.tsv", sep="\t")

print("Predictions scored:", len(allep))
print("Aggregated rows (one per mutation):", len(agg))
print("Filtered candidates:", len(filt))

cols = ["Gene Name", "Mutation", "HLA Allele", "MT Epitope Seq", "WT Epitope Seq",
        "Best MT IC50 Score", "Corresponding WT IC50 Score",
        "Corresponding Fold Change", "Tumor DNA VAF", "Tumor RNA VAF", "Gene Expression"]
top = filt[cols].sort_values("Best MT IC50 Score")
print(top.to_string(index=False))
top.to_csv("/pvacseq_output_data/top_candidates.csv", index=False)

print("\nGenes:", sorted(filt["Gene Name"].unique()))
print("Candidates per HLA allele:")
print(filt["HLA Allele"].value_counts().to_string())

# simple plot of predicted binding, if matplotlib is available
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = top["Gene Name"] + " " + top["MT Epitope Seq"]
    plt.figure(figsize=(8, 5))
    plt.barh(labels[::-1], top["Best MT IC50 Score"][::-1])
    plt.axvline(500, linestyle="--")
    plt.xlabel("Best predicted MT IC50 (nM), lower = stronger binding")
    plt.title("Filtered neoantigen candidates, HCC1395")
    plt.tight_layout()
    plt.savefig("/pvacseq_output_data/candidates_ic50.png", dpi=150)
    print("Saved plot")
except Exception as e:
    print("Plot skipped:", e)
