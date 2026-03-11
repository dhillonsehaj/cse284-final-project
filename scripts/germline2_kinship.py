import pandas as pd

# Load pairwise kinship estimates
df = pd.read_csv("results/GERMLINE2/pairwise_ibd012.csv")

# Classify relationships based on kinship coefficient
df["relationship"] = "Unrelated"
df.loc[df["kinship_phi_from_P"] > 0.044, "relationship"] = "3rd degree"
df.loc[df["kinship_phi_from_P"] > 0.088, "relationship"] = "2nd degree"
df.loc[df["kinship_phi_from_P"] > 0.177, "relationship"] = "1st degree"

# Print counts
print("\nRelationship counts:\n")
print(df["relationship"].value_counts())

# Save full table with classification
df.to_csv("results/pairwise_ibd012_classified.csv", index=False)

print("\nSaved classified pairs to:")
print("results/GERMLINE2/pairwise_ibd012_classified.csv")