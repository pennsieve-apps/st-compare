# Third-party software

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following third-party
components, installed from PyPI at the pinned versions in `requirements.txt`; their licences
apply to the image. Nothing here is GPL (the image has no scanpy, igraph or leidenalg), so the
image as a whole can be redistributed under these permissive licences.

| Component | What | Licence | Source |
|---|---|---|---|
| PyDESeq2 0.5.4 | DESeq2's negative binomial GLM, dispersion and fold-change estimation, Wald test, Cook's and independent filtering, in Python. | MIT | https://github.com/owkin/PyDESeq2 |
| formulaic 1.2.2, formulaic-contrasts 1.0.0 | Design matrices and contrasts for PyDESeq2. | MIT / BSD-3-Clause | https://github.com/matthewwardrop/formulaic, https://github.com/scverse/formulaic-contrasts |
| anndata 0.13.4 | Reading the `obs` / `var` index of `.h5ad` files (`anndata.io.read_elem`). | BSD-3-Clause | https://github.com/scverse/anndata |
| h5py 3.16.0 | Reading the count matrix of `.h5ad` files a block of rows at a time. | BSD-3-Clause | https://www.h5py.org |
| NumPy 2.5.4 | Arrays. | BSD-3-Clause (+ bundled 0BSD, MIT, Zlib, CC0-1.0 parts) | https://numpy.org |
| SciPy 1.18.1 | Sparse matrices; Welch t, Mann-Whitney U, ANOVA, Kruskal-Wallis (`scipy.stats`). | BSD-3-Clause | https://scipy.org |
| pandas 3.1.0 | Tables. | BSD-3-Clause | https://pandas.pydata.org |
| scikit-learn 1.9.1 | Dependency of PyDESeq2. | BSD-3-Clause | https://scikit-learn.org |
| Matplotlib 3.11.2 | PNG figures. | Matplotlib licence (PSF-based, BSD-compatible) | https://matplotlib.org |
| array-api-compat, cloudpickle, contourpy, cycler, donfig, fonttools, google-crc32c, interface_meta, joblib, kiwisolver, legacy-api-wrap, msgspec, narwhals, natsort, numcodecs, packaging, pillow, pydantic (+ pydantic-core, pydantic-settings, annotated-types, typing-inspection), pyparsing, python-dateutil, python-dotenv, PyYAML, scverse-misc, session-info, session-info2, six, stdlib-list, threadpoolctl, typing_extensions, tzdata, wrapt, zarr | Dependencies of the above (versions in `requirements.txt`). | MIT / BSD-2/3-Clause / MIT-CMU (Pillow) / Apache-2.0 (google-crc32c, tzdata) / Apache-2.0 or BSD-2-Clause (packaging, python-dateutil) / MPL-2.0 (legacy-api-wrap, session-info2; unmodified) / PSF-2.0 (typing_extensions) | PyPI |
| awslambdaric 2.2.1 (+ simplejson 4.2.0) | AWS Lambda runtime interface client, for the Lambda runtime. | Apache-2.0 (simplejson: MIT / AFL-2.1) | https://github.com/aws/aws-lambda-python-runtime-interface-client |
| python:3.12.16-slim-bookworm | Base image (Debian 12 + CPython 3.12). | PSF (Python); Debian packages under their own free licences | https://hub.docker.com/_/python |

## Not bundled

- **DESeq2** (R/Bioconductor, LGPL-3.0) is the method PyDESeq2 reimplements; no R code is in
  the image.
- **scCODA** (https://github.com/theislab/scCODA, BSD-3-Clause) is the fuller Bayesian model
  for compositional changes mentioned in the README; it is not bundled or called.
- **Demo data** (not in the image or the repository): Dawo S, Nonchev K, Silina K. 10x Visium
  Spatial Transcriptomics Dataset: Kidney (3) and Lung (5) Cancer with Tertiary Lymphoid
  Structures. Zenodo 2025. doi:10.5281/zenodo.14620362, CC BY 4.0.
