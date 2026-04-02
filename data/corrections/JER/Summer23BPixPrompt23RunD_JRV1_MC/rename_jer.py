from pathlib import Path

folder = Path(".")

for f in folder.glob("*.txt"):

    name = f.name

    # Skip already renamed
    if name.endswith(".jr.txt") or name.endswith(".jersf.txt"):
        continue

    if "_SF_" in name:
        new_name = f.with_suffix(".jersf.txt")

    elif "Resolution" in name:
        new_name = f.with_suffix(".jr.txt")

    else:
        print(f"Skipping (unknown): {name}")
        continue

    print(f"{name} -> {new_name.name}")
    f.rename(new_name)

print("Done.")