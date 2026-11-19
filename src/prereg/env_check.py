"""Records the computing environment and compares it with plan A, Section 7 / preregistration/v1_materials/requirements_lock.txt."""
import os
import platform
import sys

import numpy, pandas, scipy, statsmodels

EXPECTED = {"python": "3.13.12", "numpy": "2.4.4", "pandas": "2.3.3", "scipy": "1.17.1", "statsmodels": "0.14.6"}


def environment():
    got = {"python": platform.python_version(), "numpy": numpy.__version__, "pandas": pandas.__version__,
           "scipy": scipy.__version__, "statsmodels": statsmodels.__version__}
    return dict(versions=got, expected=EXPECTED, matches_plan=got == EXPECTED, platform=platform.platform(),
                executable_name=os.path.basename(sys.executable))
