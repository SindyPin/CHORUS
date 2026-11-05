"""CHORUS: Cross cohort Harmonisation of Omics for Robust and Unified early detection Signatures.

The package implements a five stage framework for pooling independent Long COVID omics
cohorts and evaluating whether pooling improves machine learning detection:

1. ``chorus.schema``      harmonised clinical metadata schema and validation
2. ``chorus.harmonise``   per layer technical harmonisation operators with fit/transform
                          semantics that respect leave-one-study-out evaluation
3. ``chorus.features``    cohort-as-environment invariant feature selection and causal priors
4. ``chorus.integrate``   mosaic multi-omics integration (late fusion, factor models)
5. ``chorus.eval``        leave-one-study-out and nested cross validation, metrics, reports

All operators are label free with respect to the held-out study, so that a harmonisation
step never leaks the outcome of the cohort used for external validation.
"""

from .io import Cohort, CohortSet  # noqa: F401
from .schema import ClinicalSchema, validate_metadata  # noqa: F401

__version__ = "0.1.0"
__all__ = ["Cohort", "CohortSet", "ClinicalSchema", "validate_metadata", "__version__"]
