#!/usr/bin/env python3
"""Index figure and table references in the manuscript source without changing it."""
import argparse
import csv
import hashlib
import re
from pathlib import Path

ASSET = re.compile(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}")
LABEL = re.compile(r"\\label\{((?:fig|tab):[^}]+)\}")
BEGIN = re.compile(r"\\begin\{(figure\*?|table\*?)\}")
END = re.compile(r"\\end\{(figure\*?|table\*?)\}")


def uncomment(line):
    return re.split(r"(?<!\\)%", line, maxsplit=1)[0]


def candidate(path):
    name = Path(path).name
    if name.startswith("limit_gaussian_s1"):
        return "scripts/gaussian_process_s1_visualization.R"
    if "vmf" in name and "grid10x10" in name:
        return "scripts/plot_paper_vmf_convergence_from_rds.R"
    if name.startswith("cycle23_joint_spatial_window"):
        return "scripts/plot_spherical_parametric_hdr_bands.R"
    if name.startswith("risoe_"):
        return "real_data/wind/plot_risoe_extended_r_density_panels.R"
    if name.endswith("simplex_logistic_gaussian_contours.pdf"):
        return "scripts/run_sediments_simplex_contours.R"
    if name.startswith("sc_vs_fmgp"):
        return "real_data/sunspots/correlacion_FMGP_p_values.r"
    return ""


def inventory(tex):
    lines = tex.read_text().splitlines()
    mode = "main"
    stack = []
    in_comment = False
    rows = []
    for number, raw in enumerate(lines, 1):
        line = uncomment(raw)
        if "\\begin{comment}" in line:
            in_comment = True
            continue
        if "\\end{comment}" in line:
            in_comment = False
            continue
        if in_comment:
            continue
        if "\\ifsupplement" in line and line.strip() == "\\ifsupplement":
            mode = "supplement"
        if match := BEGIN.search(line):
            stack.append(dict(kind=match.group(1), mode=mode, labels=[], assets=[]))
        labels = LABEL.findall(line)
        if stack:
            stack[-1]["labels"].extend(labels)
            stack[-1]["assets"].extend((number, path) for path in ASSET.findall(line))
        else:
            for label in labels:
                if label.startswith("tab:"):
                    rows.append((mode, "table", label, number, "", "", ""))
            for path in ASSET.findall(line):
                rows.append((mode, "figure", "", number, path, str((tex.parent / path).exists()).lower(), candidate(path)))
        if END.search(line) and stack:
            block = stack.pop()
            label = "|".join(block["labels"])
            if block["assets"]:
                for asset_line, path in block["assets"]:
                    rows.append((block["mode"], "figure", label, asset_line, path,
                                 str((tex.parent / path).exists()).lower(), candidate(path)))
            elif "table" in block["kind"]:
                rows.extend((block["mode"], "table", x, number, "", "", "") for x in block["labels"])
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tex", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("paper/manuscript_inventory.tsv"))
    args = parser.parse_args()
    rows = inventory(args.tex)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("mode", "kind", "label", "tex_line", "asset", "asset_exists", "producer_candidate"))
        writer.writerows(rows)
    print(f"Indexed {len(rows)} references from SHA-256 {hashlib.sha256(args.tex.read_bytes()).hexdigest()}")
