"""Per layer technical harmonisation operators.

Every operator follows scikit-learn style ``fit(X, studies, y=None)`` / ``transform(X, studies)``
semantics with one extra rule: ``fit`` only ever sees training studies, and ``transform`` on
a study that was not seen at fit time must work using label free information from that
study alone (the "addon" setting of Hornung et al. 2017). This is what makes a
leave-one-study-out estimate honest. Operators that cannot do this raise.
"""

from .operators import (  # noqa: F401
    HARMONISERS,
    AddonComBat,
    FeatureQuantileMap,
    GenePairRank,
    Identity,
    QuantileToReference,
    StudyCenterScale,
    WithinSampleRank,
    make_harmoniser,
)
from .anchored import InvarianceAnchored, sign_agreement  # noqa: F401

HARMONISERS["chorus"] = InvarianceAnchored
# chorus with the neutrophil and lymphocyte marker scores protected in the joint loading fit
HARMONISERS["chorus_protected"] = lambda **kw: InvarianceAnchored(**{"protect_scores": ["neutrophil", "lymphocyte"], **kw})

from .transcriptomics import filter_low_counts, log_cpm, select_tissue  # noqa: F401
from .proteomics import olink_overlap, within_study_zscore, ms_median_normalise  # noqa: F401
