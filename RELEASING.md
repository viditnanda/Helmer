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
3. The workflow builds the source distribution, checks that it builds and passes the tests on Ubuntu, attaches
   it to a GitHub release, and uploads it to PyPI.

Only a source distribution is published: `pip install helmer-whitney` compiles it, so FLINT and msolve must
be installed (README, *Install*). Binary wheels would have to bundle FLINT, GMP, MPFR and msolve, built from
source in the manylinux and macOS images; that is the next step for distribution.
