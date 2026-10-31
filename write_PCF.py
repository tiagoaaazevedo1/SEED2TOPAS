"""Write the OpenTOPAS parameter control files of one implant.

  - TsDicomPatient: synthetic 1 mm image on the TPS dose grid + anonymised copy of the RT Structure Set. Each structure
    gets its material with MaterialByRTStructNames (first match wins: urethra, prostate, rectum).
  - Seeds: one parallel world over the patient (Layered Mass Geometry). Geant4 ignores layered materials unless another
    parallel world with a scorer exists (tested with Ti and Pb capsules), hence CheckGrid/GridCheck.
  - Dose: unfiltered DoseToMedium in the patient. TOPAS's scoring filters lose 40-75% of the dose with layered seeds
    (tested), so the structure dose is Dose inside the Mask_<structure> voxels.
  - Every event carries one history of each source, so Mean = Sum / max(N_i), and with N_i proportional to S_K_i:
    D(Gy) = Mean x max(N_i) x 0.01 sum(S_K) / (N_total x S_K per history) x tau.
Media follow TG-186 for ultrasound plans with no CT (Beaulieu et al. 2012, sec. IV.B.2.a and Table III).
"""
import datetime
import json
import os
import shutil

import numpy as np
import pydicom
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

from plan_data import canonical

HERE = os.path.dirname(os.path.abspath(__file__))

TISSUES = {
    # Woodard & White 1986, as listed in TG-186 Table III
    "ProstateTissue": ({"H": 10.5, "C": 8.9, "N": 2.5, "O": 77.4, "Na": 0.2, "P": 0.1, "S": 0.2, "K": 0.2}, 1.04),
    # ICRU Report 46, as listed in TG-186 Table III
    "SoftTissue": ({"H": 10.5, "C": 25.6, "N": 2.7, "O": 60.2, "Na": 0.1, "P": 0.2, "S": 0.3, "Cl": 0.2, "K": 0.2}, 1.03),
    "Water": ({"H": 11.1894, "O": 88.8106}, 0.998),
}
MEDIA = {  # structure -> composition, first match wins; "Outside" is everything else
    "tissue": {"Urethra": "SoftTissue", "Prostate": "ProstateTissue", "Rectum": "SoftTissue", "Outside": "SoftTissue"},
    "water": {"Urethra": "Water", "Prostate": "Water", "Rectum": "Water", "Outside": "Water"},
}
G4 = {"H": "G4_H", "C": "G4_C", "N": "G4_N", "O": "G4_O", "Na": "G4_Na", "P": "G4_P", "S": "G4_S", "Cl": "G4_Cl", "K": "G4_K"}
PHYSICS = ['sv:Ph/Default/Modules = 3 "g4em-livermore" "g4decay" "g4radioactivedecay"',
           "d:Ph/Default/SetProductionCutLowerEdge = 1 keV", 'b:Ph/Default/Fluorescence = "True"', 'b:Ph/Default/Auger = "True"']


def scoring_grid(imp, margin=10.0, step=1.0):
    """1 mm grid over structures and seeds + margin, with voxel centres on the TPS dose grid."""
    pts = np.vstack([p for n, p in imp.structures.items() if n in ("Prostate", "Urethra", "Rectum", "PTV")] + [imp.seeds])
    lo, hi = pts.min(0) - margin, pts.max(0) + margin
    lo = imp.grid_ref + np.floor((lo - imp.grid_ref) / step) * step
    hi = imp.grid_ref + np.ceil((hi - imp.grid_ref) / step) * step
    axes = [np.arange(a, b + step / 2, step) for a, b in zip(lo, hi)]
    return dict(zip("xyz", axes))


def material_block(name, comp, rho):
    els, tot = list(comp), sum(comp.values())
    return (f'b:Ma/{name}/BuildFromMaterials = "True"\n'
            f'sv:Ma/{name}/Components = {len(els)} ' + " ".join(f'"{G4[e]}"' for e in els) + "\n"
            f'uv:Ma/{name}/Fractions = {len(els)} ' + " ".join(f"{comp[e] / tot:.6f}" for e in els) + "\n"
            f"d:Ma/{name}/Density = {rho} g/cm3\n")


def _anonymise(ds, label):
    def strip(d, elem):
        if elem.VR == "PN":
            del d[elem.tag]
    ds.walk(strip)
    for tag in ("PatientBirthDate", "PatientSex", "PatientAge", "OtherPatientIDs", "InstitutionName", "InstitutionAddress",
                "AccessionNumber", "StudyID", "StudyDescription", "SeriesDescription"):
        if tag in ds:
            delattr(ds, tag)
    ds.PatientName, ds.PatientID = label, label
    today = datetime.date.today().strftime("%Y%m%d")
    for tag in ("StudyDate", "SeriesDate", "StructureSetDate", "InstanceCreationDate"):
        if tag in ds:
            setattr(ds, tag, today)
    return ds


def write_dicom(rs_path, grid, outdir, label):
    """Synthetic 1 mm image (all voxels 0) on the scoring grid + anonymised RT Structure Set with English ROI names."""
    os.makedirs(outdir, exist_ok=True)
    ss = pydicom.dcmread(rs_path)
    for_uid = (ss.ReferencedFrameOfReferenceSequence[0].FrameOfReferenceUID
               if "ReferencedFrameOfReferenceSequence" in ss else ss.FrameOfReferenceUID)
    ss.FrameOfReferenceUID = for_uid
    study, series = generate_uid(), generate_uid()
    ss.StudyInstanceUID = study
    for r in ss.StructureSetROISequence:
        r.ROIName = canonical(r.ROIName)
    _anonymise(ss, label).save_as(os.path.join(outdir, "RS_structures.dcm"))
    nx, ny = len(grid["x"]), len(grid["y"])
    today = datetime.date.today().strftime("%Y%m%d")
    for k, z in enumerate(grid["z"]):
        meta = FileMetaDataset()
        meta.MediaStorageSOPClassUID = pydicom.uid.CTImageStorage
        meta.MediaStorageSOPInstanceUID = generate_uid()
        meta.TransferSyntaxUID = ExplicitVRLittleEndian
        ds = Dataset()
        ds.file_meta = meta
        ds.SOPClassUID, ds.SOPInstanceUID = meta.MediaStorageSOPClassUID, meta.MediaStorageSOPInstanceUID
        ds.Modality, ds.PatientName, ds.PatientID = "CT", label, label
        ds.StudyInstanceUID, ds.SeriesInstanceUID, ds.FrameOfReferenceUID = study, series, for_uid
        ds.StudyDate = ds.SeriesDate = today
        ds.InstanceNumber = k + 1
        ds.ImagePositionPatient = [float(grid["x"][0]), float(grid["y"][0]), float(z)]
        ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
        ds.PixelSpacing = [1.0, 1.0]
        ds.SliceThickness = 1.0
        ds.Rows, ds.Columns = ny, nx
        ds.SamplesPerPixel, ds.PhotometricInterpretation = 1, "MONOCHROME2"
        ds.BitsAllocated, ds.BitsStored, ds.HighBit, ds.PixelRepresentation = 16, 16, 15, 1
        ds.RescaleIntercept, ds.RescaleSlope = 0, 1
        ds.PixelData = np.zeros((ny, nx), np.int16).tobytes()
        ds.save_as(os.path.join(outdir, f"CT{k + 1:03d}.dcm"), enforce_file_format=True)


def write_implant(imp, rs_path, model, medium, outdir, label, n_total, threads, rng_seed=1, margin=10.0):
    os.makedirs(outdir, exist_ok=True)
    grid = scoring_grid(imp, margin)
    write_dicom(rs_path, grid, os.path.join(outdir, "dicom"), label)
    media = MEDIA[medium]
    mat = lambda s: f"{s}_medium"
    centre = np.array([(grid[k][0] + grid[k][-1]) / 2 for k in "xyz"])
    n_i = np.round(n_total * imp.sk / imp.sk.sum()).astype(int)
    structs = [s for s in media if s != "Outside" and s in imp.structures]
    sk_h = model.sk_per_hist
    gy = None if not sk_h else float(n_i.max() * 0.01 * imp.sk.sum() / (n_i.sum() * sk_h) * imp.tau_h)
    L = [f"# SEED2TOPAS: {len(imp.seeds)} x {model.name}, {medium}",
         "# Dose.csv is Gy per event: " + (f"total dose (Gy) = value x {gy:.6e}" if gy else
                                           "no S_K calibration for this seed yet, see ../seed_calibration/"),
         f"i:Ts/NumberOfThreads = {threads}", f"i:Ts/Seed = {rng_seed}",
         f"i:Ts/ShowHistoryCountAtInterval = {max(1, int(n_i.max()) // 20)}",
         'b:Ts/PauseBeforeQuit = "False"', "i:Ts/ParameterizationErrorMaxReports = 1"] + PHYSICS + [
         'sv:Ph/Default/LayeredMassGeometryWorlds = 1 "Seeds"', "",
         model.materials]
    L += [material_block(mat(s), *TISSUES[comp]) for s, comp in media.items()]
    L += ["d:Ge/World/HLX = 30 cm", "d:Ge/World/HLY = 30 cm", "d:Ge/World/HLZ = 30 cm", 's:Ge/World/Material = "Vacuum"',
          "# surrounding medium beyond the patient image; the image is centred in it",
          's:Ge/Body/Type = "TsBox"', 's:Ge/Body/Parent = "World"', f's:Ge/Body/Material = "{mat("Outside")}"',
          "d:Ge/Body/HLX = 10 cm", "d:Ge/Body/HLY = 10 cm", "d:Ge/Body/HLZ = 10 cm",
          "# ---- patient: 1 mm image on the dose grid + the RT Structure Set",
          's:Ge/Patient/Type = "TsDicomPatient"', 's:Ge/Patient/Parent = "Body"', 's:Ge/Patient/DicomDirectory = "dicom"',
          's:Ge/Patient/ImagingToMaterialConverter = "MaterialTagNumber"',
          "iv:Ge/Patient/MaterialTagNumbers = 1 0", f'sv:Ge/Patient/MaterialNames = 1 "{mat("Outside")}"',
          "# material of each structure (first match wins where structures overlap)",
          f"sv:Ge/Patient/MaterialByRTStructNames = {len(structs)} " + " ".join(f'"{s}"' for s in structs),
          f"sv:Ge/Patient/MaterialByRTStructMaterials = {len(structs)} " + " ".join(f'"{mat(s)}"' for s in structs),
          "# ---- seeds at the plan positions, in one parallel world over the patient (Layered Mass Geometry)",
          's:Ge/Seeds/Type = "Group"', 's:Ge/Seeds/Parent = "Patient"', 'b:Ge/Seeds/IsParallel = "True"']
    for i, (s, n) in enumerate(zip(imp.seeds, n_i), start=1):
        x, y, z = s - centre
        L.append(model.block(i, x, y, z, n, "Seeds", True))
    hl = [len(grid[k]) / 2 for k in "xyz"]
    L += ["# A second parallel world with a scorer is required: without one the layered seed materials are ignored.",
          's:Ge/CheckGrid/Type = "TsBox"', 's:Ge/CheckGrid/Parent = "Patient"', 'b:Ge/CheckGrid/IsParallel = "True"',
          f"d:Ge/CheckGrid/HLX = {hl[0]} mm", f"d:Ge/CheckGrid/HLY = {hl[1]} mm", f"d:Ge/CheckGrid/HLZ = {hl[2]} mm",
          f"i:Ge/CheckGrid/XBins = {len(grid['x'])}", f"i:Ge/CheckGrid/YBins = {len(grid['y'])}", f"i:Ge/CheckGrid/ZBins = {len(grid['z'])}",
          's:Sc/GridCheck/Quantity = "DoseToMedium"', 's:Sc/GridCheck/Component = "CheckGrid"',
          's:Sc/GridCheck/OutputType = "csv"', 's:Sc/GridCheck/OutputFile = "GridCheck"',
          's:Sc/GridCheck/IfOutputFileAlreadyExists = "Overwrite"', 'sv:Sc/GridCheck/Report = 1 "Mean"',
          "# ---- dose in the whole patient",
          's:Sc/Dose/Quantity = "DoseToMedium"', 's:Sc/Dose/Component = "Patient"', 's:Sc/Dose/OutputType = "csv"',
          's:Sc/Dose/OutputFile = "Dose"', 's:Sc/Dose/IfOutputFileAlreadyExists = "Overwrite"',
          'sv:Sc/Dose/Report = 2 "Mean" "Standard_Deviation"',
          "# ---- structure masks as TOPAS reads them (-1 = outside). Use only the mask: filtered dose is not reliable here."]
    scored = [s for s in ("Prostate", "Urethra", "Rectum") if s in imp.structures]
    for s in scored:
        L += [f's:Sc/Mask_{s}/Quantity = "DoseToMedium"', f's:Sc/Mask_{s}/Component = "Patient"',
              f'sv:Sc/Mask_{s}/OnlyIncludeIfInRTStructure = 1 "{s}"',
              f'b:Sc/Mask_{s}/SetBinToMinusOneIfNotInRTStructure = "True"',
              f's:Sc/Mask_{s}/OutputType = "csv"', f's:Sc/Mask_{s}/OutputFile = "Mask_{s}"',
              f's:Sc/Mask_{s}/IfOutputFileAlreadyExists = "Overwrite"', f'sv:Sc/Mask_{s}/Report = 1 "Sum"']
    open(os.path.join(outdir, "implant.txt"), "w").write("\n".join(L) + "\n")
    json.dump({"seeds_mm": imp.seeds.tolist(), "sk_U": imp.sk.tolist(), "half_life_d": imp.half_life_d, "tau_h": imp.tau_h,
               "prescription_gy": imp.prescription_gy, "source_name": imp.source_name, "seed_model": model.name, "seed_file": model.file,
               "sk_per_hist_Gy_cm2": sk_h, "scoring_grid": {k: v.tolist() for k, v in grid.items()},
               "centre_mm": centre.tolist(), "histories_per_seed": n_i.tolist(), "n_total": int(n_i.sum()), "medium": medium,
               "media": {s: [mat(s), c] for s, c in media.items()}, "scored_structures": scored, "gy_per_event": gy},
              open(os.path.join(outdir, "implant.json"), "w"))
    write_view(outdir, model, len(imp.seeds), scored)
    return gy


def write_view(outdir, model, nseeds, scored):
    """view.txt: seeds only, no histories, scorers off. Parallel worlds and patient voxels would draw as white grids."""
    vols = [f"{c}_{i}" for i in range(1, nseeds + 1) for c in model.visible]
    L = ["# Visualisation only (SEED2TOPAS): no histories, scorers off. Run: topas view.txt", "includeFile = implant.txt",
         "i:Ts/NumberOfThreads = 1", 'b:Ts/UseQt = "True"', 'b:Ts/PauseBeforeQuit = "True"',
         f"sv:Gr/ViewA/VisibleWorlds = {len(vols)} " + " ".join(f'"{v}"' for v in vols),
         's:Gr/ViewA/Type = "OpenGL"', "i:Gr/ViewA/WindowSizeX = 1000", "i:Gr/ViewA/WindowSizeY = 900",
         "d:Gr/ViewA/Theta = 60 deg", "d:Gr/ViewA/Phi = 30 deg", 'b:Gr/ViewA/IncludeAxes = "True"', "d:Gr/ViewA/AxesSize = 2 cm",
         'b:Ge/QuitIfOverlapDetected = "False"', 'b:Ge/World/Invisible = "True"', 'b:Ge/Body/Invisible = "True"',
         'b:Ge/CheckGrid/Invisible = "True"', 'b:Ge/Patient/Invisible = "True"']
    L += [f"ic:So/{model.source}_{i}/NumberOfHistoriesInRun = 0" for i in range(1, nseeds + 1)]
    L += [f'b:Sc/{s}/Active = "False"' for s in ["Dose", "GridCheck"] + [f"Mask_{x}" for x in scored]]
    open(os.path.join(outdir, "view.txt"), "w").write("\n".join(L) + "\n")


def write_calibration(model, outdir, n=1_000_000_000):
    """Air-kerma strength per history of a single seed in vacuum: WAFAC-like detector (7.6 deg half-angle) at 100 cm,
    track-length estimator with air mu_en/rho, photons above 5 keV. S_K per history = Mean x (100 cm)^2."""
    os.makedirs(outdir, exist_ok=True)
    shutil.copy(os.path.join(HERE, "data", "Muen.dat"), outdir)
    L = [f"# SEED2TOPAS: air-kerma strength of {model.name} (WAFAC geometry, in vacuum). After the run:",
         "#   S_K per history (Gy cm2) = Mean in kermaLTE.csv x 100 x 100; add to the seed file: # sk_per_hist_Gy_cm2 = <value>",
         "i:Ts/NumberOfThreads = 8", "i:Ts/Seed = 1", f"i:Ts/ShowHistoryCountAtInterval = {max(1, n // 20)}",
         'b:Ts/PauseBeforeQuit = "False"'] + PHYSICS + ["d:Ph/Default/CutForElectron = 10 cm",
         "d:Ge/World/HLX = 110 cm", "d:Ge/World/HLY = 110 cm", "d:Ge/World/HLZ = 110 cm", 's:Ge/World/Material = "Vacuum"',
         model.materials,
         model.block(1, 0.0, 0.0, 0.0, n, "World", False),
         's:Ge/WAFAC/Type = "TsSphere"', 's:Ge/WAFAC/Parent = "World"', 's:Ge/WAFAC/Material = "G4_AIR"',
         "d:Ge/WAFAC/RMin = 99.95 cm", "d:Ge/WAFAC/RMax = 100.05 cm", "d:Ge/WAFAC/SPhi = 0 deg", "d:Ge/WAFAC/DPhi = 360 deg",
         "d:Ge/WAFAC/STheta = 82.4 deg", "d:Ge/WAFAC/DTheta = 15.2 deg",
         's:Sc/Kerma/Quantity = "TrackLengthEstimator"', 's:Sc/Kerma/Component = "WAFAC"', 's:Sc/Kerma/InputFile = "Muen.dat"',
         's:Sc/Kerma/OutputType = "csv"', 's:Sc/Kerma/OutputFile = "kermaLTE"', 's:Sc/Kerma/IfOutputFileAlreadyExists = "Overwrite"',
         'sv:Sc/Kerma/Report = 2 "Mean" "Standard_Deviation"', 'sv:Sc/Kerma/OnlyIncludeParticlesNamed = 1 "gamma"',
         "d:Sc/Kerma/OnlyIncludeIfIncidentParticleKEAbove = 5 keV"]
    open(os.path.join(outdir, "kerma_WAFAC.txt"), "w").write("\n".join(L) + "\n")
