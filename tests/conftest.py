import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
for sub in ("workflow/lib", "workflow/scripts"):
    sys.path.insert(0, str(REPO / sub))
