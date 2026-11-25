# Releasing Helmer

`.github/workflows/release.yml` does the work once two one-time steps are done.

## One-time setup

1. **PyPI.** The name `helmer` is taken on PyPI by an unrelated project, so the distribution is
   `helmer-whitney` (the module is still `import helmer`). At https://pypi.org, under *Your account →
   Publishing → Add a new pending publisher*, enter: project `helmer-whitney`, owner `viditnanda`,
   repository `Helmer`, workflow `release.yml`, environment `pypi`.
2. **Zenodo (DOI).** At https://zenodo.org (logged in with GitHub), under *GitHub*, switch on
   `viditnanda/Helmer`. Every GitHub release then gets a DOI; `.zenodo.json` and `CITATION.cff` supply the
   metadata.

## Each release

1. Update the version in `pyproject.toml`, `python/helmer/__init__.py` and `CITATION.cff` (and its
   `date-released`), and date the top entry of `CHANGELOG.md`.
2. Tag and push: `git tag v0.1.0 && git push origin v0.1.0`.
3. The workflow builds the binary wheels (`wheels.yml`: Linux x86-64 and macOS on Apple silicon and Intel,
   CPython 3.9-3.13, each tested in a clean environment) and the source distribution, attaches them all to a
   GitHub release, and uploads them to PyPI.

The wheels bundle GMP, MPFR, FLINT and msolve, built by `tools/build_wheel_deps.sh` without CPU-specific
instructions (see the comments there), so `pip install helmer-whitney` needs nothing else on those
platforms. Elsewhere (Windows, Linux on ARM) pip builds the source distribution, which needs FLINT and
msolve installed. `wheels.yml` also runs on every push to a branch named `wheels`, and by hand.
