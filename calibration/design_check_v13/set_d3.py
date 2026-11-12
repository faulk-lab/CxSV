#!/usr/bin/env python3
"""Write the D3 depth-noise thresholds into config/scoring.yaml (text edit; comments kept)."""
import re, sys
p = sys.argv[1]
s = open(p).read()
old_block = re.search(r"# Thresholds: PLACEHOLDERS.*?default_sigma_max: null[^\n]*\n", s, re.S).group(0)
new_block = """# Thresholds (panel-v1.3 design check D3, docs/design_check_v13/): per PoN stratum,
# median + 3.5 x scaled MAD (robust z 3.5) of the CN sigma of the 1,619 panel-v1.2-r2
# libraries (cn/stats.json). Panel libraries above: ONT_1kGP_1019 22/1,017 (mostly
# < 5x depth), PacBio_HiFi 10/135, ONT_LRSC_R10 0/247, ONT_LRSC_R9 0/220. Nine test
# genomes (all R10): BOC_5883 (0.322) and BOC_7308 (0.237) above, the others 0.16-0.21.
# Every PoN stratum is listed, so default_sigma_max (strata without a PoN) stays null.
depth_noise:
  enabled: true
  sigma_max:
    ONT_LRSC_R10: 0.233    # median 0.1629, scaled MAD 0.0200, n 247
    ONT_LRSC_R9: 0.220     # median 0.1558, scaled MAD 0.0185, n 220
    ONT_1kGP_1019: 0.401   # median 0.2330, scaled MAD 0.0480, n 1,017
    PacBio_HiFi: 0.198     # median 0.1586, scaled MAD 0.0111, n 135
  default_sigma_max: null     # strata not listed above
"""
s = s.replace(old_block, new_block)
open(p, "w").write(s)
import yaml
print(yaml.safe_load(open(p))["depth_noise"])
