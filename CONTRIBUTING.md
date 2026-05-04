# Contributing to SED_Model

Thank you for your interest in contributing to SED_Model.

SED_Model is a scientific software package for forward synthetic photometry and inverse stellar-parameter inference from atmosphere grids and photometric filters.

## Reporting Bugs & Requesting Features

Open an issue on GitHub. For bugs, include:

- A clear description of the problem
- Steps to reproduce it
- What you expected vs. what happened
- Operating system
- Python version
- Compiler version, if the issue involves the Fortran extension
- Relevant grid/filter names or paths, if applicable
- Full traceback or test failure output

For feature requests, describe the scientific use case and why it belongs in SED_Model rather than SED_Tools.

## Contribution Workflow

1. Fork the repository.
2. Create a branch for your change.
3. Make your changes.
4. Run the relevant tests.
5. Open a pull request against `main` with a clear description of what changed and why.

Keep PRs focused: one logical change per PR.

Large API changes, changes to the Fortran kernels, or changes to photometric zero-point conventions should be discussed in an issue before opening a PR.

## Development Setup

```bash
git clone https://github.com/nialljmiller/SED_Model.git
cd SED_Model
python -m pip install -e ".[dev]"
python setup.py build_ext --inplace
````

Run tests with:

```bash
pytest tests/test_forward.py
pytest tests/test_inverse.py
```

## Code Style

There is no strict formatter yet. Match the style of the surrounding code.

Prefer:

* Clear scientific variable names
* Explicit units in comments/docstrings
* Small functions with testable behaviour
* Reproducible examples
* Minimal hidden global state

## What to Contribute

Useful contributions include:

* Bug fixes with a clear root cause
* Additional tests for forward/inverse consistency
* New extinction-law validation
* Performance improvements to interpolation or filter convolution
* Better documentation and examples
* Robust handling of additional SED_Tools grid products

## What Not to Do

* Do not bundle unrelated changes into one PR.
* Do not submit numerical changes without tests or a clear reference.
* Do not change zero-point definitions without discussion.
* Do not change the public API casually; this package is intended to be citeable and stable.
* Do not commit generated files such as `*.so`, `*.pyc`, `build/`, `dist/`, `*.egg-info/`, large grids, or downloaded data.

## Relationship to SED_Tools

SED_Model is downstream of SED_Tools. SED_Tools manages SED/filter data; SED_Model performs modelling and inference using those data.

Changes that affect data acquisition, standardization, or grid construction likely belong in SED_Tools. Changes that affect synthetic magnitudes, forward modelling, extinction, likelihoods, or posterior inference likely belong in SED_Model.

## Contact

For questions not suited to a public issue: **[niall.j.miller@gmail.com](mailto:niall.j.miller@gmail.com)**
