"""SEED2TOPAS: parametrize a VariSeed seed implant (RT Structure Set, RT Dose, RT Plan) into OpenTOPAS parameter control
files, with the seeds of the plan and the structures inside the simulation. NOT for clinical use.

    python3 SEED2TOPAS.py --m gui
    python3 SEED2TOPAS.py --m inputfile input.txt
"""
import os
import shutil
import sys

import plan_data
import seed_model
import write_PCF as W

MEDIA = ("tissue", "water", "both")
MODELS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "seed models")


def print_help():
    print(__doc__)
    print("input file, one value per line:")
    for k, t in enumerate(("Project name (no patient name or ID: it goes into the output)", "path to DICOM-structure file (RS)",
                           "path to DICOM-dose file (RD)", "path to DICOM-RT plan file (RP)",
                           'path to the seed geometry file (a bare file name is looked up in "seed models")',
                           'Medium: "tissue", "water" or "both"', "Total histories, e.g. 1e9", "Threads", "Output folder"), start=1):
        print(f"  ({k}) {t}")


def seed_path(p):
    p = os.path.expanduser(p.strip())
    if not os.path.isfile(p) and os.path.isfile(os.path.join(MODELS, p)):
        return os.path.join(MODELS, p)
    return p


def convert(project, rs, rd, rp, seed_file, medium, histories, threads, out, log=print):
    """Writes <out>/<project>/<project>_<medium>/ and run.sh. Returns the project folder."""
    project = project.strip()
    if not project or os.sep in project:
        raise ValueError("invalid project name")
    if medium not in MEDIA:
        raise ValueError(f"medium must be one of {MEDIA}")
    seed_file = seed_path(seed_file)
    for name, p in (("RS", rs), ("RD", rd), ("RP", rp), ("Seed geometry", seed_file)):
        if not os.path.isfile(p):
            raise FileNotFoundError(f"{name} file not found: {p}")
    imp = plan_data.load(rs, rd, rp)
    if not {"Prostate", "Urethra", "Rectum"} & set(imp.structures):
        raise ValueError(f"no prostate, urethra or rectum in the RS (ROIs: {', '.join(sorted(imp.structures))})")
    low = project.lower()
    for i in imp.identifiers:
        for token in str(i).replace("^", " ").split():
            if len(token) > 2 and token.lower() in low:
                raise ValueError("the project name contains the patient's name or ID: choose another one")
    n_total = int(float(histories))
    if n_total < len(imp.seeds):
        raise ValueError(f"at least {len(imp.seeds)} histories are needed (one per seed)")

    model = seed_model.SeedModel(seed_file)
    warnings = list(imp.warnings) + model.notes
    if model.matches(imp.source_name) is False:
        warnings.append(f'the RP says "{imp.source_name}" but the seed file is {model.name}')
    if model.half_life_d and abs(imp.half_life_d - model.half_life_d) > 0.5:
        warnings.append(f"half-life in the RP {imp.half_life_d} d, seed file {model.half_life_d} d. The RP value is used")
    if model.status != "validated":
        warnings.append(f"{model.file} is not marked as validated: check its TG-43 parameters before trusting doses")

    root = os.path.join(out, project)
    if os.path.exists(os.path.join(root, "run.sh")):
        raise FileExistsError(f"{root} already has a generated project: choose another name or move it")
    media = ("tissue", "water") if medium == "both" else (medium,)
    folders = []
    log(f"{len(imp.seeds)} seeds, {imp.source_name}, S_K {imp.sk.mean():.3f} U, total {imp.sk.sum():.2f} U, "
        f"T1/2 {imp.half_life_d} d, Rx {imp.prescription_gy:g} Gy")
    log(f"Structures: {', '.join(sorted(imp.structures))}")
    log(f"Seed: {model.name} [{model.status}], from {model.file}")
    for med in media:
        name = f"{project}_{med}"
        gy = W.write_implant(imp, rs, model, med, os.path.join(root, name), project, n_total, int(threads))
        folders.append(name)
        log(f"  {name}: Gy per event x " + (f"{gy:.6e}" if gy else "(needs the S_K calibration)"))
    shutil.copy(seed_file, os.path.join(root, "seed_model_used.txt"))
    lines = ["#!/usr/bin/env bash", "# Each run overwrites its outputs. If topas is not on the PATH: TOPAS=/path/to/topas bash run.sh",
             'cd "$(dirname "$0")"', 'TOPAS="${TOPAS:-topas}"']
    lines += [f'(cd "{f}" && "$TOPAS" implant.txt)' for f in folders]
    if not model.sk_per_hist:
        W.write_calibration(model, os.path.join(root, "seed_calibration"))
        lines.insert(4, '(cd seed_calibration && "$TOPAS" kerma_WAFAC.txt)')
        warnings.append(f"{model.file} has no S_K per history, so doses cannot be converted to Gy: run seed_calibration/kerma_WAFAC.txt, "
                        "then add \"# sk_per_hist_Gy_cm2 = <Mean x 100 x 100>\" to the seed file and regenerate")
    with open(os.path.join(root, "run.sh"), "w") as f:
        f.write("\n".join(lines) + "\n")
    os.chmod(os.path.join(root, "run.sh"), 0o755)
    for w in warnings:
        log("WARNING: " + w)
    log(f"\nWritten to {root}\nRun:  bash \"{os.path.join(root, 'run.sh')}\"\nLook first:  cd \"{os.path.join(root, folders[0])}\" && topas view.txt")
    return root


def from_file(path):
    v = [ln.strip() for ln in open(path) if ln.strip()]
    if len(v) != 9:
        raise ValueError(f"{path} must have 9 lines, found {len(v)}")
    convert(*v)


def gui():
    import tkinter as tk
    from tkinter import filedialog, ttk
    root = tk.Tk()
    root.title("SEED2TOPAS")
    frm = ttk.Frame(root, padding=12)
    frm.grid()
    vals = {}

    def row(r, label, key, default="", browse=None):
        ttk.Label(frm, text=label).grid(row=r, column=0, sticky="w", pady=2)
        vals[key] = tk.StringVar(value=default)
        ttk.Entry(frm, textvariable=vals[key], width=60).grid(row=r, column=1, pady=2)
        if browse:
            ttk.Button(frm, text="...", width=3, command=lambda: vals[key].set(browse() or vals[key].get())).grid(row=r, column=2)

    def pick():
        return filedialog.askopenfilename(filetypes=[("DICOM", "*.dcm"), ("All files", "*")])

    def pick_seed():
        return filedialog.askopenfilename(initialdir=MODELS, title="Seed geometry file", filetypes=[("TOPAS", "*.txt"), ("All files", "*")])

    models = sorted(f for f in os.listdir(MODELS) if f.endswith(".txt")) if os.path.isdir(MODELS) else []
    row(0, "Project name (no patient data)", "project", "case01")
    row(1, "RT Structure Set (RS)", "rs", browse=pick)
    row(2, "RT Dose (RD)", "rd", browse=pick)
    row(3, "RT Plan (RP)", "rp", browse=pick)
    row(4, "Seed geometry file", "seed", os.path.join(MODELS, models[0]) if len(models) == 1 else "", browse=pick_seed)
    ttk.Label(frm, text="Medium").grid(row=5, column=0, sticky="w", pady=2)
    vals["medium"] = tk.StringVar(value="both")
    ttk.Combobox(frm, textvariable=vals["medium"], values=MEDIA, state="readonly").grid(row=5, column=1, sticky="w")
    row(6, "Total histories", "hist", "1e9")
    row(7, "Threads", "threads", "8")
    row(8, "Output folder", "out", os.path.expanduser("~"), browse=filedialog.askdirectory)
    text = tk.Text(frm, width=90, height=14)
    text.grid(row=10, column=0, columnspan=3, pady=8)

    def log(msg):
        text.insert("end", msg + "\n")
        text.see("end")

    def run():
        text.delete("1.0", "end")
        try:
            convert(vals["project"].get(), vals["rs"].get(), vals["rd"].get(), vals["rp"].get(), vals["seed"].get(),
                    vals["medium"].get(), vals["hist"].get(), vals["threads"].get(), vals["out"].get(), log)
        except Exception as e:
            log(f"ERROR: {e}")

    ttk.Button(frm, text="Write TOPAS files", command=run).grid(row=9, column=1, pady=6)
    root.mainloop()


def main():
    a = sys.argv[1:]
    if a == ["--m", "gui"]:
        gui()
    elif len(a) == 3 and a[:2] == ["--m", "inputfile"]:
        from_file(a[2])
    else:
        print_help()
        sys.exit(0 if a in ([], ["--help"]) else 1)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError, FileExistsError) as e:
        sys.exit(f"ERROR: {e}")
