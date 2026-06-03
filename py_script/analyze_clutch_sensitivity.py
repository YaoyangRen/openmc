#!/usr/bin/env python3
"""Analyze OpenMC generalized CLUTCH sensitivity HDF5 output.

Default usage:

    python analyze_clutch_sensitivity.py

Optional:

    python analyze_clutch_sensitivity.py clutch_sensitivity.h5 --csv result.csv
    python analyze_clutch_sensitivity.py clutch_sensitivity.h5 --tree
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np


@dataclass
class Parameter:
    pid: int
    variable: str
    material_id: int
    nuclide: str
    value: float


@dataclass
class MethodResult:
    name: str
    available: bool
    dlogk: np.ndarray
    sensitivity: np.ndarray
    uncertainty: np.ndarray
    numerator: np.ndarray
    denominator: float
    batch_ids: np.ndarray
    batch_numerator: np.ndarray
    batch_denominator: np.ndarray
    uses_transfer_function: bool = False


def decode_value(value: Any) -> Any:
    """Decode HDF5 scalar/string values into ordinary Python objects."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.bytes_):
        return bytes(value).decode("utf-8")
    if isinstance(value, np.ndarray):
        if value.dtype.kind in {"S", "O", "U"}:
            return [decode_value(x) for x in value.tolist()]
        return value
    return value


def read_string_dataset(group: h5py.Group, name: str) -> list[str]:
    if name not in group:
        return []
    data = group[name][...]
    return [str(decode_value(x)) for x in data]


def read_scalar_attr(attrs: h5py.AttributeManager, name: str, default: Any = None) -> Any:
    if name not in attrs:
        return default
    return decode_value(attrs[name])


def read_parameters(h5: h5py.File) -> list[Parameter]:
    if "parameters" not in h5:
        raise KeyError("Missing /parameters group")

    params = h5["parameters"]
    ids = params["ids"][...].astype(int)
    variables = read_string_dataset(params, "variable")
    material_ids = params["material_id"][...].astype(int)
    nuclides = read_string_dataset(params, "nuclide")
    values = params["parameter_value"][...].astype(float)

    result = []
    for i, pid in enumerate(ids):
        result.append(
            Parameter(
                pid=int(pid),
                variable=variables[i] if i < len(variables) else "",
                material_id=int(material_ids[i]),
                nuclide=nuclides[i] if i < len(nuclides) else "",
                value=float(values[i]),
            )
        )
    return result


def _read_optional_array(group: h5py.Group, name: str, default_shape: tuple[int, ...]) -> np.ndarray:
    if name in group:
        return np.asarray(group[name][...])
    return np.zeros(default_shape)


def read_method(group: h5py.Group, name: str, n_params: int) -> MethodResult:
    method = group[name]
    available = bool(int(read_scalar_attr(method.attrs, "available", 1)))
    denominator_data = np.asarray(method["denominator"][...], dtype=float)
    denominator = float(denominator_data.flat[0]) if denominator_data.size else 0.0

    dlogk = _read_optional_array(method, "dlogk_dparameter", (n_params,)).astype(float)
    sensitivity = _read_optional_array(method, "sensitivity", (n_params,)).astype(float)
    numerator = _read_optional_array(method, "numerator", (n_params,)).astype(float)
    uncertainty = _read_optional_array(method, "uncertainty", (n_params,)).astype(float)
    batch_ids = _read_optional_array(method, "batch_ids", (0,)).astype(int)
    batch_denominator = _read_optional_array(method, "batch_denominator", (0,)).astype(float)
    batch_numerator = _read_optional_array(
        method, "batch_numerator", (batch_denominator.size, n_params)
    ).astype(float)

    return MethodResult(
        name=name,
        available=available,
        dlogk=dlogk,
        sensitivity=sensitivity,
        uncertainty=uncertainty,
        numerator=numerator,
        denominator=denominator,
        batch_ids=batch_ids,
        batch_numerator=batch_numerator,
        batch_denominator=batch_denominator,
        uses_transfer_function=bool(int(read_scalar_attr(method.attrs, "uses_transfer_function", 0))),
    )


def read_methods(h5: h5py.File, n_params: int) -> dict[str, MethodResult]:
    if "method" not in h5:
        raise KeyError("Missing /method group")
    methods = {}
    method_group = h5["method"]
    for name in sorted(method_group.keys()):
        if isinstance(method_group[name], h5py.Group):
            methods[name] = read_method(method_group, name, n_params)
    return methods


def print_header(h5: h5py.File, path: Path) -> None:
    print("=" * 96)
    print("CLUTCH sensitivity analysis")
    print("=" * 96)
    print(f"file: {path}")
    for key in (
        "filetype",
        "version",
        "sensitivity_type",
        "primary_method",
        "method",
        "adjoint_source",
        "cclutch_method",
    ):
        value = read_scalar_attr(h5.attrs, key)
        if value is not None:
            print(f"{key}: {value}")
    print()


def print_diagnostics(h5: h5py.File) -> None:
    if "diagnostics" not in h5:
        return

    diag = h5["diagnostics"]
    print("Diagnostics")
    print("-" * 96)
    for key in (
        "total_fission_sites",
        "total_scored_sites",
        "total_dropped_sites",
        "total_cclutch_events",
        "total_cclutch_scored_events",
        "total_cclutch_dropped_events",
        "total_cclutch_missing_source_events",
        "active_batches_scored",
    ):
        value = read_scalar_attr(diag.attrs, key)
        if value is not None:
            print(f"{key}: {value}")

    total = read_scalar_attr(diag.attrs, "total_fission_sites", 0)
    scored = read_scalar_attr(diag.attrs, "total_scored_sites", 0)
    dropped = read_scalar_attr(diag.attrs, "total_dropped_sites", 0)
    if total:
        print(f"scored fraction: {float(scored) / float(total):.6f}")
        print(f"dropped fraction: {float(dropped) / float(total):.6f}")

    for name in ("grid_shape", "grid_lower_left", "grid_upper_right", "grid_pitch"):
        if name in diag:
            print(f"{name}: {np.asarray(diag[name][...]).tolist()}")
    print()


def format_rel_unc(value: float, unc: float) -> str:
    if value == 0.0:
        return "nan"
    return f"{abs(unc / value):.3e}"


def print_method_table(params: list[Parameter], method: MethodResult) -> None:
    label = method.name
    if method.name == "cclutch_history" and method.uses_transfer_function:
        label += " (transfer-function C-CLUTCH)"
    print(f"Method: {label}")
    print("-" * 120)
    print(
        f"{'id':>5} {'variable':>16} {'mat':>6} {'nuclide':>12} "
        f"{'parameter':>14} {'dlogk/dp':>14} {'sensitivity':>14} "
        f"{'unc':>12} {'rel_unc':>12}"
    )
    print("-" * 120)
    for i, param in enumerate(params):
        print(
            f"{param.pid:5d} {param.variable:>16} {param.material_id:6d} "
            f"{param.nuclide or '-':>12} {param.value:14.6e} "
            f"{method.dlogk[i]:14.6e} {method.sensitivity[i]:14.6e} "
            f"{method.uncertainty[i]:12.4e} "
            f"{format_rel_unc(method.sensitivity[i], method.uncertainty[i]):>12}"
        )
    print(f"denominator mean: {method.denominator:.8e}")
    print()


def print_batch_diagnostics(method: MethodResult) -> None:
    den = method.batch_denominator
    if den.size == 0:
        print(f"Batch diagnostics for {method.name}: no batch data")
        print()
        return

    print(f"Batch diagnostics: {method.name}")
    print("-" * 96)
    valid = den != 0.0
    zero_count = int(np.count_nonzero(~valid))
    print(f"batches: {den.size}")
    print(f"zero denominator batches: {zero_count}")
    print(f"denominator min/mean/max: {den.min():.6e} / {den.mean():.6e} / {den.max():.6e}")
    if den.size > 1:
        print(f"denominator std: {den.std(ddof=1):.6e}")
        if den.mean() != 0.0:
            print(f"denominator coeff. variation: {den.std(ddof=1) / abs(den.mean()):.6e}")

    if method.batch_numerator.ndim == 2 and method.batch_numerator.size and np.any(valid):
        ratios = np.full_like(method.batch_numerator, np.nan, dtype=float)
        ratios[valid, :] = method.batch_numerator[valid, :] / den[valid, None]
        ratio_mean = np.nanmean(ratios, axis=0)
        ratio_std = np.nanstd(ratios, axis=0, ddof=1) if np.count_nonzero(valid) > 1 else np.zeros(ratios.shape[1])
        print("batch ratio dlogk/dp mean/std by parameter:")
        for i, (mean, std) in enumerate(zip(ratio_mean, ratio_std), start=1):
            print(f"  param[{i:02d}]: mean={mean:.6e}, std={std:.6e}")
    print()


def print_method_comparison(params: list[Parameter], methods: dict[str, MethodResult]) -> None:
    if "fclutch_fm" not in methods or "cclutch_history" not in methods:
        return

    primary = methods["fclutch_fm"]
    diag = methods["cclutch_history"]
    print("Method comparison: C-CLUTCH - F-CLUTCH")
    print("-" * 96)
    print(f"{'id':>5} {'variable':>16} {'delta_sens':>14} {'rel_delta':>14}")
    print("-" * 96)
    for i, param in enumerate(params):
        delta = diag.sensitivity[i] - primary.sensitivity[i]
        rel = delta / primary.sensitivity[i] if primary.sensitivity[i] != 0.0 else np.nan
        print(f"{param.pid:5d} {param.variable:>16} {delta:14.6e} {rel:14.6e}")
    print()


def print_h5_tree(group: h5py.Group, prefix: str = "/") -> None:
    for key in sorted(group.keys()):
        item = group[key]
        path = f"{prefix}{key}"
        if isinstance(item, h5py.Dataset):
            print(f"{path}: dataset shape={item.shape}, dtype={item.dtype}")
        elif isinstance(item, h5py.Group):
            print(f"{path}/: group")
            print_h5_tree(item, f"{path}/")


def write_csv(path: Path, params: list[Parameter], methods: dict[str, MethodResult]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "method",
                "id",
                "variable",
                "material_id",
                "nuclide",
                "parameter_value",
                "dlogk_dparameter",
                "sensitivity",
                "uncertainty",
                "relative_uncertainty",
                "numerator",
                "denominator",
            ]
        )
        for method in methods.values():
            for i, param in enumerate(params):
                rel_unc = (
                    abs(method.uncertainty[i] / method.sensitivity[i])
                    if method.sensitivity[i] != 0.0
                    else np.nan
                )
                writer.writerow(
                    [
                        method.name,
                        param.pid,
                        param.variable,
                        param.material_id,
                        param.nuclide,
                        param.value,
                        method.dlogk[i],
                        method.sensitivity[i],
                        method.uncertainty[i],
                        rel_unc,
                        method.numerator[i],
                        method.denominator,
                    ]
                )


def analyze(filename: Path, csv_path: Path | None = None, show_tree: bool = False) -> None:
    if not filename.exists():
        raise FileNotFoundError(f"Cannot find {filename}")

    with h5py.File(filename, "r") as h5:
        print_header(h5, filename)
        if show_tree:
            print("HDF5 tree")
            print("-" * 96)
            print_h5_tree(h5)
            print()

        params = read_parameters(h5)
        methods = read_methods(h5, len(params))

        print_diagnostics(h5)
        for method_name in ("fclutch_fm", "cclutch_history"):
            if method_name in methods:
                print_method_table(params, methods[method_name])
                print_batch_diagnostics(methods[method_name])

        for method_name, method in methods.items():
            if method_name not in {"fclutch_fm", "cclutch_history"}:
                print_method_table(params, method)
                print_batch_diagnostics(method)

        print_method_comparison(params, methods)

        if csv_path is not None:
            write_csv(csv_path, params, methods)
            print(f"CSV written: {csv_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze OpenMC clutch_sensitivity.h5 output."
    )
    parser.add_argument(
        "filename",
        nargs="?",
        default="clutch_sensitivity.h5",
        help="Path to clutch_sensitivity.h5",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Optional CSV output path for method/parameter table.",
    )
    parser.add_argument(
        "--tree",
        action="store_true",
        help="Print the HDF5 tree before the analysis tables.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    analyze(Path(args.filename), csv_path=args.csv, show_tree=args.tree)


if __name__ == "__main__":
    main()
