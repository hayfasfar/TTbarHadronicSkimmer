from pathlib import Path

# change this if needed
folder = Path(".")

for f in folder.glob("*.txt"):

    name = f.name

    # Skip already renamed files
    if name.endswith(".jec.txt") or name.endswith(".junc.txt"):
        continue

    # -------- uncertainty files --------
    if any(key in name for key in [
        "Uncertainty",
        "UncertaintySources",
        "Regrouped"
    ]):
        new_name = f.with_suffix(".junc.txt")

    # -------- correction files --------
    else:
        new_name = f.with_suffix(".jec.txt")

    print(f"{f.name}  ->  {new_name.name}")

    f.rename(new_name)

print("Done.")