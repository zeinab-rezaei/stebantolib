

import argparse
from pathlib import Path
from typing import List, Tuple

# splash import comes from the spectra‑hash package
from splash import Splash
from splash import Spectrum
from splash import SpectrumType
from typing import Callable, Any  # new import at the top

from collections import Counter
import re
import json, time, os ,argparse
import requests
from urllib.parse import quote
from pathlib import Path




def splash_for_peaks(peaks: List[Tuple[float, float]]) -> str:
    """Generate a SPLASH identifier from a list of (mz, intensity)."""
    # Spectra‑hash expects intensities scaled to arbitrary units – raw is fine.
    if len(peaks) == 0 :
        raise ValueError(f"Record starting with has no peaks.")
        print(f"Record starting with {current.get('COMPOUND_NAME','<unknown>')} has no peaks.")
    spectrum = Spectrum(peaks, SpectrumType.MS)
    return Splash().splash(spectrum)

