# This repository is in experimental stage!

This is intended to maybe become the next generation course material
for teaching Dune/PDELab. For now, it only documents some of my
experiments with jupyter/xeus-cling/binder etc.

Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dokempf/dune-jupyter-course/master)

## Development

To join the development make sure to read this section.

### Pre-commit

Before contributing, you should run this:

```
python -m pip install pre-commit
pre-commit install
```

It will make sure that the notebook outputs do not end up in the repository erroneously.

### Extensions

The notebooks use the unofficial jupyter notebook extensions, an installation guide can be found [here](https://jupyter-contrib-nbextensions.readthedocs.io/en/latest/install.html).
Currently, the following extensions are in use:
 - Collapsible Headings
 - Equation Auto Numbering
 - Exercise2
 - Table of Contents (2)
 - Rubberband (required by Exercise2)
