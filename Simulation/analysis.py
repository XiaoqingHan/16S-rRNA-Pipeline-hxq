import pandas as pd, glob, os

rows = []
for run in glob.glob("runs/*/summary.csv"):
    df = pd.read_csv(run)
    df["run"] = os.path.basename(os.path.dirname(run))
    rows.append(df)

if not rows:
    print("[implementationanalysis] No runs/*/summary.csv found. Nothing to aggregate.")
else:
    all = pd.concat(rows, ignore_index=True)
    all.to_csv("all_runs_summary.csv", index=False)
    print(f"[implementationanalysis] Wrote all_runs_summary.csv with {len(all)} rows.")