# conda-forge (not a feedstock yet)

The in-repo recipe is [`conda/meta.yaml`](meta.yaml). It installs the
**Python** package (`flexaidds`), not the native engine.

This is **not** `conda install -c conda-forge flexaidds`. That channel does
not carry the package until a feedstock is accepted. First-shot today:

```bash
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git
cd FlexAIDdS
conda env create -f environment.yml
conda activate flexaidds
python -c "import flexaidds as fd; print(fd.__version__)"
```

To publish on conda-forge later: open a staged-recipes PR using this
`meta.yaml` (version still derived from `python/flexaidds/__version__.py`).
Do not document `conda-forge` as a working first-shot until that feedstock
exists and `conda search -c conda-forge flexaidds` returns a version.
