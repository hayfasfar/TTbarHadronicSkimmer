import zipfile
import sys
from pathlib import Path

zip_path = Path(sys.argv[1])
out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else zip_path.stem

out_dir.mkdir(parents=True, exist_ok=True)

with zipfile.ZipFile(zip_path, 'r') as z:
    z.extractall(out_dir)

print(f"Extracted to {out_dir}")