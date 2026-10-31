"""Read the VariSeed RT Plan (seeds, S_K), RT Dose (grid) and RT Structure Set (contours), and check they belong together."""
import dataclasses
import math
import unicodedata

import numpy as np
import pydicom

ALIASES = {"Prostate": ("prostate", "prostata", "prost", "ctv"), "Urethra": ("urethra", "uretra", "uret"),
           "Rectum": ("rectum", "reto", "recto", "rect"), "PTV": ("ptv",)}


def canonical(name):
    """English ASCII structure name, as TOPAS needs it (Próstata -> Prostate)."""
    n = unicodedata.normalize("NFKD", str(name)).encode("ascii", "ignore").decode().strip().lower()
    for canon, keys in ALIASES.items():
        if n in keys:
            return canon
    return str(name).strip()


@dataclasses.dataclass
class Implant:
    seeds: np.ndarray            # n x 3, mm, DICOM patient coordinates
    sk: np.ndarray               # U, per seed, at the reference date of the plan
    half_life_d: float
    source_name: str             # as written by the TPS
    prescription_gy: float
    grid_ref: np.ndarray         # first voxel centre of the TPS dose grid, mm
    structures: dict             # {canonical name: (points n x 3, mm)}
    identifiers: list            # patient name and ID, only to check they do not reach the output
    warnings: list

    @property
    def tau_h(self):
        """Mean life, h: total dose of a permanent implant = initial dose rate x tau."""
        return self.half_life_d * 24.0 / math.log(2.0)


def _read(path, modality):
    ds = pydicom.dcmread(path)
    if ds.Modality != modality:
        raise ValueError(f"{path}\n  is {ds.Modality}, expected {modality}")
    return ds


def _frame(ds):
    if "ReferencedFrameOfReferenceSequence" in ds:
        return ds.ReferencedFrameOfReferenceSequence[0].FrameOfReferenceUID
    return getattr(ds, "FrameOfReferenceUID", None)


def load(rs_path, rd_path, rp_path):
    warn = []
    plan, dose, ss = _read(rp_path, "RTPLAN"), _read(rd_path, "RTDOSE"), _read(rs_path, "RTSTRUCT")

    src = {int(s.SourceNumber): s for s in plan.SourceSequence}
    seeds, sk = [], []
    for app in plan.ApplicationSetupSequence:
        for ch in app.ChannelSequence:
            cps = ch.BrachyControlPointSequence
            p0 = np.array(cps[0].ControlPoint3DPosition, float)
            if not np.allclose(p0, np.array(cps[-1].ControlPoint3DPosition, float)):
                raise ValueError("the plan has a moving source: only permanent seed implants are supported")
            seeds.append(p0)
            sk.append(float(src[int(ch.ReferencedSourceNumber)].ReferenceAirKermaRate))
    if not seeds:
        raise ValueError("no seeds in the RT Plan")
    s0 = plan.SourceSequence[0]
    rx = [float(r.TargetPrescriptionDose) for r in getattr(plan, "DoseReferenceSequence", [])
          if hasattr(r, "TargetPrescriptionDose")]

    if not np.allclose(np.array(dose.ImageOrientationPatient, float), [1, 0, 0, 0, 1, 0]):
        raise ValueError("the RT Dose grid is not axial (identity orientation expected)")

    structures = {}
    names = {r.ROINumber: canonical(r.ROIName) for r in ss.StructureSetROISequence}
    for rc in ss.ROIContourSequence:
        pts = [np.array(c.ContourData, float).reshape(-1, 3) for c in getattr(rc, "ContourSequence", [])]
        if pts:
            structures[names[rc.ReferencedROINumber]] = np.vstack(pts)
    for need in ("Prostate", "Urethra", "Rectum"):
        if need not in structures:
            warn.append(f"no {need} contour found in the RT Structure Set (ROIs: {', '.join(sorted(structures))})")

    uids = {"RP": _frame(plan), "RD": _frame(dose), "RS": _frame(ss)}
    if len({u for u in uids.values() if u}) > 1:
        warn.append("RP, RD and RS have different Frame of Reference UIDs: check that they are the same case")

    seeds = np.array(seeds)
    for i in range(len(seeds)):       # seeds that overlap would overlap in TOPAS too
        d = seeds[i + 1:] - seeds[i]
        close = (np.hypot(d[:, 0], d[:, 1]) < 0.8) & (np.abs(d[:, 2]) < 4.5)
        for j in np.where(close)[0]:
            warn.append(f"seeds {i + 1} and {i + j + 2} overlap (seeds along z, 0.8 x 4.5 mm)")

    ids = [str(x) for x in (plan.PatientName, plan.PatientID) if str(x).strip()]
    return Implant(seeds, np.array(sk), float(s0.SourceIsotopeHalfLife), str(s0.SourceIsotopeName), rx[0] if rx else float("nan"),
                   np.array(dose.ImagePositionPatient, float), structures, ids, warn)
