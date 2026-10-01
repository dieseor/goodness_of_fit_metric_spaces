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


def table_rows(table):
    """Split a tabular block into rows of stripped cells, ignoring comments and rules."""
    body = "\n".join(re.split(r"(?<!\\)%", line, maxsplit=1)[0] for line in table.splitlines())
    rows = []
    for row in re.split(r"\\\\", body):
        row = re.sub(r"\\(toprule|midrule|bottomrule|hline)\b|\\cmidrule(\([^)]*\))?\{[^}]*\}", "", row)
        if "&" in row:
            rows.append([cell.strip().split("\n")[-1].strip() for cell in row.split("&")])
    return rows


def label(cell):
    """Plain text of a row label such as \\texttt{WhiteCells\\_microscopic} or $\\mathrm{C}_2$."""
    text = re.sub(r"\\(texttt|mathrm|textbf|mathbf)\{([^}]*)\}", r"\2", cell)
    text = text.replace("\\_", "_")
    text = re.sub(r"_\{?(\d)\}?", r"\1", text)
    return re.sub(r"[${}\\\s]", "", text).replace("--", "_").lower()


def header_index(rows, names):
    for i, row in enumerate(rows):
        if all(name in row for name in names):
            return i
    raise AssertionError(f"Header with {names} not found")


def compare(where, table_value, saved, digits):
    expected = f"{float(saved):.{digits}f}"
    found = number(table_value)
    if found is None or f"{found:.{digits}f}" != expected:
        raise AssertionError(f"{where}: table has {table_value!r}, saved result is {expected}")


def check_real_data(tables):
    """Compare every saved value with its own row and column of the real-data tables."""
    checked = 0

    # Sunspots: a single row of estimates followed by the KS and CM p-values.
    t = table_rows(tables["sunspots"])
    h = header_index(t, ["KS", "CM"])
    data = [row for row in t[h + 1:] if len(row) == len(t[h])]
    assert len(data) == 1, data
    saved = {row["statistic_type"].split("_")[0]: row for row in rows("sunspots_gof.csv")}
    estimates = ["temporal_weight1", "temporal_alpha1", "temporal_beta1", "temporal_alpha2", "temporal_beta2",
                 "a_N", "b_N", "c"]
    for position, column in enumerate(estimates):
        compare(f"sunspots {column}", data[0][position], saved["ks"][column], 2)
        checked += 1
    for column, statistic in (("KS", "ks"), ("CM", "cvm")):
        compare(f"sunspots {column}", data[0][t[h].index(column)], saved[statistic]["p_value"], 4)
        checked += 1

    # Comets: models in rows, long and short period in column blocks.
    t = table_rows(tables["comets"])
    h = header_index(t, ["AIC", "BIC", "KS", "CM"])
    periods = ("long", "short") if " ".join(t[h - 1]).find("Long") < " ".join(t[h - 1]).find("Short") else ("short", "long")
    ks_cols = [i for i, cell in enumerate(t[h]) if cell == "KS"]
    cm_cols = [i for i, cell in enumerate(t[h]) if cell == "CM"]
    saved = {(row["model"].lower(), row["period"]): (row["ks_pvalue"], row["cvm_pvalue"]) for row in rows("comets_c2_sc.csv")}
    saved.update({("ub", row["dataset"].split("_")[0]): (row["gof_ks_p_value"], row["gof_cvm_p_value"])
                  for row in rows("comets_ub.csv")})
    body = [row for row in t[h + 1:] if len(row) == len(t[h])]
    assert sorted(label(row[0]) for row in body) == ["c2", "sc", "ub"], [row[0] for row in body]
    for row in body:
        for period, ks_col, cm_col in zip(periods, ks_cols, cm_cols):
            ks, cm = saved[(label(row[0]), period)]
            compare(f"comets {row[0]} {period} KS", row[ks_col], ks, 4)
            compare(f"comets {row[0]} {period} CM", row[cm_col], cm, 4)
            checked += 2

    # Wind: one row per window.
    t = table_rows(tables["wind"])
    h = header_index(t, ["Months", "KS", "CM"])
    saved = {row["window_id"]: row for row in rows("wind_gof.csv")}
    body = [row for row in t[h + 1:] if len(row) == len(t[h])]
    assert sorted(label(row[0]) for row in body) == sorted(saved), [row[0] for row in body]
    for row in body:
        s = saved[label(row[0])]
        compare(f"wind {row[0]} n", row[1], s["n_valid"], 0)
        mu = [float(x) for x in re.findall(r"-?\d+\.\d+", row[2])]
        assert [f"{x:.2f}" for x in mu] == [f"{float(s[k]):.2f}" for k in ("mu1_hat", "mu2_hat", "mu3_hat")], (row[0], mu)
        compare(f"wind {row[0]} kappa", row[3], s["kappa_hat"], 1)
        compare(f"wind {row[0]} KS", row[t[h].index("KS")], s["p_value_KS"], 4)
        compare(f"wind {row[0]} CM", row[t[h].index("CM")], s["p_value_CvM"], 4)
        checked += 7

    # Compositions: one row per dataset, with the BHEP column.
    t = table_rows(tables["compositions"])
    h = header_index(t, ["Dataset", "KS", "CM", "BHEP"])
    saved = {row["dataset"].lower(): row for row in rows("compositions_ks_cvm.csv")}
    bhep = {row["dataset"].lower(): row["p_value"] for row in rows("hz_pvalues.csv")}
    body = [row for row in t[h + 1:] if len(row) == len(t[h])]
    assert sorted(label(row[0]) for row in body) == sorted(saved) == sorted(bhep), [row[0] for row in body]
    for row in body:
        s = saved[label(row[0])]
        compare(f"compositions {row[0]} n", row[1], s["n"], 0)
        compare(f"compositions {row[0]} d", row[2], s["ilr_dimension"], 0)
        compare(f"compositions {row[0]} KS", row[t[h].index("KS")], s["ks_pvalue"], 4)
        compare(f"compositions {row[0]} CM", row[t[h].index("CM")], s["cvm_pvalue"], 4)
        compare(f"compositions {row[0]} BHEP", row[t[h].index("BHEP")], bhep[label(row[0])], 4)
        checked += 5
    return checked


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

    cells = check_real_data(tables)
    print(f"Matched 480 simulation cells, all 80 Table 3 values and {cells} cells of the four real-data tables.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tex", type=Path)
    main(parser.parse_args().tex)
