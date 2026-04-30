from __future__ import annotations

import argparse
import copy
import glob
from pathlib import Path
from typing import Any, Iterable

from coffea import processor, util
from hist import Hist


def _as_paths(paths: Iterable[str | Path]) -> list[Path]:
    selected = [Path(path).expanduser() for path in paths]
    missing = [str(path) for path in selected if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing coffea files:\n" + "\n".join(missing))
    if not selected:
        raise ValueError("No input coffea files were provided.")
    return selected


def _add_mapping_values(left: dict[Any, Any], right: dict[Any, Any]) -> dict[Any, Any]:
    result = copy.deepcopy(left)
    for key, value in right.items():
        if key in result:
            result[key] = result[key] + value
        else:
            result[key] = copy.deepcopy(value)
    return result


def _combine_value(key: str, left: Any, right: Any, source: Path) -> Any:
    if isinstance(left, Hist) and isinstance(right, Hist):
        try:
            return left + right
        except Exception as exc:
            raise ValueError(f"Histogram {key!r} is not compatible in {source}") from exc

    if isinstance(left, processor.defaultdict_accumulator) and isinstance(
        right, processor.defaultdict_accumulator
    ):
        return _add_mapping_values(left, right)

    if isinstance(left, processor.dict_accumulator) and isinstance(
        right, processor.dict_accumulator
    ):
        return _add_mapping_values(left, right)

    if isinstance(left, processor.list_accumulator) and isinstance(
        right, processor.list_accumulator
    ):
        return left + right

    return left


def combine_coffea_outputs(
    input_files: Iterable[str | Path],
    output_file: str | Path | None = None,
) -> dict[str, Any]:
    """Combine top-level histograms and accumulator counters from coffea outputs.

    The first file is used as the template. Compatible ``hist.Hist`` objects are
    added bin-by-bin, cutflow/weight-style accumulators are added by key, and
    metadata-like objects such as ``analysisCategories`` are kept from the first
    file. If ``output_file`` is provided, the combined object is written with
    ``coffea.util.save``.
    """

    paths = _as_paths(input_files)
    combined = copy.deepcopy(util.load(paths[0]))
    if not isinstance(combined, dict):
        raise TypeError(f"{paths[0]} did not load to a dictionary-like coffea output.")

    for source in paths[1:]:
        current = util.load(source)
        if not isinstance(current, dict):
            raise TypeError(f"{source} did not load to a dictionary-like coffea output.")

        for key, value in current.items():
            if key not in combined:
                combined[key] = copy.deepcopy(value)
                continue
            combined[key] = _combine_value(key, combined[key], value, source)

    combined["combined_inputs"] = [str(path) for path in paths]

    if output_file is not None:
        output_path = Path(output_file).expanduser()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        util.save(combined, output_path)

    return combined


def _expand_input_patterns(patterns: Iterable[str]) -> list[Path]:
    files: list[Path] = []
    for pattern in patterns:
        if any(ch in pattern for ch in "*?[]"):
            matches = sorted(Path(path) for path in glob.glob(pattern))
            files.extend(matches)
        else:
            files.append(Path(pattern))
    return files


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Combine QCD coffea outputs by summing their histograms and counters."
    )
    parser.add_argument("inputs", nargs="+", help="Input .coffea files or glob patterns")
    parser.add_argument("-o", "--output", required=True, help="Output .coffea file")
    args = parser.parse_args()

    input_files = _expand_input_patterns(args.inputs)
    combined = combine_coffea_outputs(input_files, args.output)
    hist_count = sum(isinstance(value, Hist) for value in combined.values())
    print(f"Wrote {args.output}")
    print(f"Combined {len(input_files)} files and {hist_count} histograms")


if __name__ == "__main__":
    main()
