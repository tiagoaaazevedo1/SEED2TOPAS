# SEED2TOPAS

SEED2TOPAS converts a permanent seed implant exported from VariSeed (RT Structure Set, RT Dose and RT Plan) into OpenTOPAS parameter control files. Each seed of the plan and the prostate, urethra and rectum are modeled in the simulation.

It is not for clinical use.

The tool follows the approach of the TPS2TOPAS interface, which does the same for TrueBeam plans. SEED2TOPAS does not share code with it. Please cite the original work:

Ortiz R, Sawkey D, Faddegon B, D-Kondo N, Ramos-Méndez J. An interface tool to parametrize treatment plans for the TrueBeam radiotherapy system into TOPAS parameter control files for Monte Carlo simulation. Physica Medica. 2024;124:104485. https://doi.org/10.1016/j.ejmp.2024.104485

## Requirements

- Python 3 with numpy and pydicom 3 or later (`pip install -r requirements.txt`)
- tkinter, only for the window
- [OpenTOPAS](https://github.com/OpenTOPAS), to run the output

## Usage

```
python3 SEED2TOPAS.py --m gui
python3 SEED2TOPAS.py --m inputfile input.txt
```

The input file has one value per line (see `inputfile_template.txt`):

1. Project name (it goes into the output, so do not use the patient's name or ID)
2. RS file
3. RD file
4. RP file
5. Seed geometry file (a file name alone is looked up in `seed models`)
6. Medium: `tissue`, `water` or `both`
7. Total histories, split over the seeds in proportion to their air-kerma strength
8. Threads
9. Output folder

With 1e9 histories we obtained about 2% Type A uncertainty per 1 mm voxel, in 4 to 5 h on 8 threads.

## Output

```
<project>/
  run.sh                 runs all the simulations
  seed_model_used.txt    copy of the seed file
  <project>_tissue/      one folder per medium
    implant.txt          TOPAS input
    view.txt             draws the seeds, without histories
    implant.json         seeds, air-kerma strength, grid, histories and dose factor
    dicom/               synthetic 1 mm image and anonymised RS, read by TOPAS
  seed_calibration/      only for a seed without calibration
```

To run one simulation: `cd <project>_tissue && topas implant.txt`. If `topas` is not on the PATH: `TOPAS=/path/to/topas bash run.sh`.

`Dose.csv` is in Gy per event. To get the total dose in Gy, multiply by `gy_per_event` from `implant.json`. `Mask_<structure>.csv` gives the voxels TOPAS assigns to each structure (-1 outside), and the dose to a structure is `Dose.csv` in those voxels. Voxels that contain a seed hold the dose to the seed materials, not to tissue.

## Model

The patient is a `TsDicomPatient`, made of a synthetic 1 mm image on the TPS dose grid and the RS. The RD is only used to align this grid. Where structures overlap, the urethra takes priority over the prostate and the prostate over the rectum.

The `tissue` option follows the TG-186 recommendations for plans without CT images (Beaulieu et al. 2012): the prostate has the composition of Woodard and White (1986), with a density of 1.04 g/cm³, and the urethra, rectum and surrounding tissues have the mean male soft tissue composition of ICRU Report 46, with a density of 1.03 g/cm³. The `water` option uses water with a density of 0.998 g/cm³ for all of them.

The seeds are placed at the positions in the RP, with their axis along z, in a parallel world over the patient (Layered Mass Geometry). Dose is scored with an unfiltered `DoseToMedium` scorer. We found two TOPAS behaviours with this geometry, which explain the input:

- The layered materials are ignored unless another parallel world has a scorer. This is the purpose of `CheckGrid`.
- Scoring filters lost 40 to 75% of the dose with layered seeds. This is why the structure dose comes from the masks.

The total dose of the permanent implant is

D = Mean × max(N_i) × 0.01 × ΣS_K / (N_total × S_K per history) × τ

where τ = T½ / ln 2.

## Seed geometry file

`seed models` contains our model of the IsoSeed I25.S17plus, validated against the TG-43 data of Pantelis et al. (2013).

A seed file holds plain TOPAS parameters for one seed at the origin with its axis along z: materials, geometry and one source. The names are free. The tool repeats the seed for each position in the plan and sets its parent, `IsParallel` and `NumberOfHistoriesInRun`. Other parameters (physics, scoring, the world and scored components, such as a WAFAC detector) are ignored with a warning, so an existing single-seed input can be used. A surrounding phantom is not recognised and has to be removed first.

Comment lines give the information TOPAS does not need:

```
# name = IsoSeed I25.S17plus (Eckert & Ziegler BEBIG)
# half_life_d = 59.4
# aliases = s17plus, s17 plus
# sk_per_hist_Gy_cm2 = 4.0269e-14
# status = validated
```

The half-life is only compared with the RP, and the RP value is used. The aliases (lower case, no dots) are compared with the source name in the RP. Without `sk_per_hist_Gy_cm2` the dose cannot be converted to Gy: the tool then writes `seed_calibration/kerma_WAFAC.txt`, and the value to add is the mean in `kermaLTE.csv` × 100 × 100. A new seed should have its TG-43 parameters checked before use.

`data/Muen.dat` is the file from the OpenTOPAS brachytherapy example.

## Patient data

The RP, RD and RS contain patient identifiers. In the output, the names and IDs are removed and the dates replaced, but the anonymised RS keeps the contours. Generated projects should not be shared without checking them, and `.gitignore` keeps them out of this repository.

## References

Beaulieu L, Carlsson Tedgren Å, Carrier JF, et al. Report of the Task Group 186 on model-based dose calculation methods in brachytherapy beyond the TG-43 formalism. Med Phys. 2012;39(10):6208-36. https://doi.org/10.1118/1.4747264

ICRU. Photon, electron, proton and neutron interaction data for body tissues. ICRU Report 46. 1992.

Pantelis E, Papagiannis P, Anagnostopoulos G, Baltas D. New 125I brachytherapy source IsoSeed I25.S17plus: Monte Carlo dosimetry simulation and comparison to sources of similar design. J Contemp Brachytherapy. 2013;5(4):240-9. https://doi.org/10.5114/jcb.2013.39631

Woodard HQ, White DR. The composition of body tissues. Br J Radiol. 1986;59(708):1209-18. https://doi.org/10.1259/0007-1285-59-708-1209
