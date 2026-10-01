#!/usr/bin/env python3
"""Compare saved numerical results with the current manuscript tables, read only."""
import argparse
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
LABELS = {
    "simulation": "tab:empirical_null_calibration",
    "timing": "tab:bootstrap_procedure_comparison",
    "sunspots": "tab:sunspots_cycle23_temporal_models",
    "comets": "tab:comets-gof",
    "wind": "tab:risoe-wind-gof",
    "compositions": "tab:logistic-gaussian-real-data",
}


def block(tex, label):
    start = tex.index(r"\label{" + label + "}")
    return tex[start:tex.index(r"\end{table}", start)]


def number(cell):
    match = re.search(r"-?\d+(?:\.\d+)?", cell)
    return float(match.group()) if match else None


def simulation_cells(table):
    result = {}
    scenario = dimension = None
    for line in table.splitlines():
        if "&" not in line or not line.rstrip().endswith(r"\\"):
            continue
        fields = [field.strip() for field in line.split("&")]
        if len(fields) != 13 or number(fields[2]) is None:
            continue
        if number(fields[0]) is not None:
            scenario = int(number(fields[0]))
        if number(fields[1]) is not None:
            dimension = int(number(fields[1]))
        values = [number(field) for field in fields[3:]]
        if scenario and dimension and all(value is not None for value in values):
            for i, value in enumerate(values):
                key = (scenario, dimension, (50, 100, 200, 400, 800)[i // 2],
                       number(fields[2]), ("ks", "cvm")[i % 2])
                result[key] = value
    return result


def rows(name):
    with (ROOT / "artifacts" / "tables" / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def check_rounded(table, values, digits=4):
    missing = [f"{float(value):.{digits}f}" for value in values
               if f"{float(value):.{digits}f}" not in table]
    if missing:
        raise AssertionError(f"Values missing from table: {missing}")


def main(tex_path):
    tex = tex_path.read_text()
    tables = {name: block(tex, label) for name, label in LABELS.items()}
    observed = simulation_cells(tables["simulation"])
    audit = json.loads((ROOT / "artifacts/tables/simulation_audit.json").read_text())
    expected = {tuple(cell["key"]): cell["paper"] for cell in audit["cells"]}
    for row in rows("scenario1_n800.csv"):
        for statistic in ("ks", "cvm"):
            expected[(1, int(row["d"]), 800, float(row["beta"]), statistic)] = float(row[statistic])
    assert len(expected) == len(observed) == 480, (len(expected), len(observed))
    differences = [(key, value, observed.get(key)) for key, value in expected.items()
                   if observed.get(key) != value]
    if differences:
        raise AssertionError(f"Simulation table differs: {differences[:5]}")

    timing_rows = []
    for line in tables["timing"].splitlines():
        fields = [field.strip() for field in line.split("&")]
        if len(fields) == 10 and number(fields[1]) in (2, 5):
            values = [number(field) for field in fields[2:]]
            if all(value is not None for value in values):
                timing_rows.append(values)
    assert len(timing_rows) == 10, len(timing_rows)
    slow = {(row["scenario"], int(row["d"])): row for row in rows("table3_slow_summary.csv")}
    fast = {(row["scenario"], int(row["d"])): row for row in rows("table3_fast_rates.csv")}
    for dimension_index, dimension in enumerate((2, 5)):
        for scenario_index in range(8):
            key = (f"S{scenario_index + 1}", dimension)
            for row_index, statistic in ((0, "ks"), (2, "cvm")):
                delta = float(slow[key][f"{statistic}_rejection_percent"]) - float(
                    fast[key][f"{statistic}_rejection_percent"])
                assert round(delta, 1) == timing_rows[row_index + dimension_index][scenario_index]
            fast_time = float(fast[key]["median_elapsed_seconds"])
            slow_time = float(slow[key]["median_elapsed_seconds"])
            assert round(fast_time, 1) == timing_rows[4 + dimension_index][scenario_index]
            assert round(slow_time, 1) == timing_rows[6 + dimension_index][scenario_index]
            assert round(slow_time / fast_time) == timing_rows[8 + dimension_index][scenario_index]

    check_rounded(tables["sunspots"], [row["p_value"] for row in rows("sunspots_gof.csv")])
    check_rounded(tables["comets"], [row[stat] for name in ("comets_c2_sc.csv", "comets_ub.csv")
                                         for row in rows(name)
                                         for stat in (("ks_pvalue", "cvm_pvalue") if name == "comets_c2_sc.csv"
                                                      else ("gof_ks_p_value", "gof_cvm_p_value"))])
    check_rounded(tables["wind"], [row[stat] for row in rows("wind_gof.csv")
                                       for stat in ("p_value_KS", "p_value_CvM")])
    check_rounded(tables["compositions"], [row[stat] for row in rows("compositions_ks_cvm.csv")
                                               for stat in ("ks_pvalue", "cvm_pvalue")])
    check_rounded(tables["compositions"], [row["p_value"] for row in rows("hz_pvalues.csv")])
    print("Matched 480 simulation cells, all 80 Table 3 values, and saved GOF p-values in four real-data tables, including the BHEP column.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tex", type=Path)
    main(parser.parse_args().tex)
