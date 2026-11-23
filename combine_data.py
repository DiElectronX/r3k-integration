import pandas as pd

INPUT_CSV = "../data/logs/cutflow_step1_mix.csv"
OUTPUT_CSV = "test.csv"

# Read CSV
df = pd.read_csv(INPUT_CSV)

# Select specific rows
mask_data = df["Sample Name"].str.startswith("data2022")
mask_ssdata = df["Sample Name"].str.startswith("samesigndata2022")

df_sel_data = df[mask_data]
df_sel_ssdata = df[mask_ssdata]

if df_sel_data.empty or df_sel_ssdata.empty:
    raise RuntimeError("No rows matching 'data2022' or 'ssdata2022' found")

# Columns that should NOT be summed
non_numeric_cols = ["Sample Name", "Region", "Trigger_Mode"]

# Sum numeric columns
summed_data = df_sel_data.drop(columns=non_numeric_cols).sum()
summed_ssdata = df_sel_ssdata.drop(columns=non_numeric_cols).sum()

# Build output rows
out_data = pd.DataFrame([summed_data])
out_ssdata = pd.DataFrame([summed_ssdata])

out_data.insert(0, "Trigger_Mode", df_sel_data["Trigger_Mode"].iloc[0])
out_data.insert(0, "Region", df_sel_data["Region"].iloc[0])
out_data.insert(0, "Sample Name", "data2022")

out_ssdata.insert(0, "Trigger_Mode", df_sel_ssdata["Trigger_Mode"].iloc[0])
out_ssdata.insert(0, "Region", df_sel_ssdata["Region"].iloc[0])
out_ssdata.insert(0, "Sample Name", "ssdata2022")

# Combine rows
out_df = pd.concat([out_data, out_ssdata], ignore_index=True)

# Write output once
out_df.to_csv(OUTPUT_CSV, index=False)

print(f"Wrote combined rows to {OUTPUT_CSV}")