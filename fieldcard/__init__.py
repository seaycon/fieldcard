"""FIELDCARD: explainable weekly field cards that attribute crop stress to a cause.

fieldcard.core        data-source-independent method (indices, FAO-56 water balance,
                      reference curves, cause attribution, JSON cards)
fieldcard.evaluation  season classification and the pre-registered TRL 4 criteria
fieldcard.sources     readers for open data (DWD, Thuenen crop maps, Sentinel-2 on
                      Planetary Computer); imported on demand, needs the [data] extras
"""
from .core import *  # noqa: F401,F403
from .core import __version__  # noqa: F401
from . import evaluation  # noqa: F401
