# Releasing Helmer

`.github/workflows/release.yml` does the work once two one-time steps are done.

## One-time setup

1. **PyPI.** The name `helmer` is taken on PyPI by an unrelated project, so the distribution is
   `helmer-strat` (the module is still `import helmer`). At https://pypi.org, under *Your account →
   Publishing → Add a new pending publisher*, enter: project `helmer-strat`, owner `viditnanda`,
   repository `Helmer`, workflow `release.yml`, environment `pypi`.
2. **Zenodo (DOI).** At https://zenodo.org (logged in with GitHub), under *GitHub*, switch on
   `viditnanda/Helmer`. Every GitHub release then gets a DOI; `.zenodo.json` and `CITATION.cff` supply the
   metadata.

## Each release

1. On a branch, update the version in `pyproject.toml`, `python/helmer/__init__.py`, `CITATION.cff` (and its
   `date-released`) and `web/index.html` (twice: the header link and "The current version is"), and date the
   top entry of `CHANGELOG.md`. Merge into `main`.
2. Tag the merged commit and push: `git tag vX.Y.Z && git push origin vX.Y.Z`. The tag must name the version
   in `pyproject.toml`; the workflow stops before publishing anything if it does not. Zenodo must already be
   switched on, since it archives only releases made after that.
3. The workflow builds the binary wheels (`wheels.yml`: Linux x86-64 and macOS on Apple silicon and Intel,
   CPython 3.9-3.13, each tested in a clean environment) and the source distribution, attaches them all to a
   GitHub release, and uploads them to PyPI.

The wheels bundle GMP, MPFR, FLINT and msolve, built by `tools/build_wheel_deps.sh` without CPU-specific
instructions (see the comments there), so `pip install helmer-strat` needs nothing else on those
platforms. Elsewhere (Windows, Linux on ARM) pip builds the source distribution, which needs FLINT and
msolve installed. `wheels.yml` also runs on every push to a branch named `wheels`, and by hand.
