# Validation

The reference analysis uses Python 3.12.14, matplotlib 3.10.7 and NumPy 2.5.4.
Exact inputs are recorded in `INPUT_MANIFEST.json`; the dependency closure is
recorded in `requirements-lock.txt`.

## Numerical checks

- Sixteen Python tests cover coordinate conventions, gaps, touching boundaries,
  invalid intervals, class definitions, masks, empty unions and conservation.
- Two hundred randomized track pairs are compared with an independent per-base
  reference calculation.
- One hundred deterministic 10-kb windows from the real inputs are checked
  against per-base reference maps. Windows may overlap; this is a software
  correctness check, not independent biological sampling.
- The six scientific output files reproduce byte-for-byte in a separate run.
- All seven code cells in `analysis.ipynb` execute in order in a fresh local
  Jupyter kernel, with saved figures and no error outputs.

## Explorer checks

JavaScript calculations match the Python results for all six mask/class
combinations. Render/event tests cover the four selectable mask/enhancer
combinations, all 225 matrix cells, matrix selection, coordinate validation,
empty regions, and CSV/SVG export content. These interface checks use a minimal
DOM adapter; full browser visual/responsive testing is not established by them.

The PNG figures have been visually reviewed for labels, denominators and
attribution. Hosted execution and Google Colab execution are outside the
reference validation environment.

## Repeat the checks

Install `requirements-dev.txt`, then run from the repository root:

```sh
python -m pytest -q
python -m chromatin_atlas --output validation-rerun
python scripts/check_reproduction.py validation-rerun
node tests/test_dashboard.cjs
node tests/test_interface.cjs
python scripts/execute_notebook.py --output validation-rerun/analysis.executed.ipynb
```

Scientific data products and the generated HTML are compared exactly. Notebook
timing metadata and figure pixels can vary with operating systems, fonts and
software builds. Software validation does not provide functional or clinical
validation of the regulatory annotations.
