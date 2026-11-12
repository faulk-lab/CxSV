#!/usr/bin/env python3
"""cn_params.json of a code tree, as rule cn_params_json writes it; prints the junction
collapse distance and INS length ratio (for extract_junctions.py).
  cn_params_of.py CODE OUT.json"""
import json, sys
import yaml
code, out = sys.argv[1:3]
y = yaml.safe_load(open(f"{code}/config/calling.yaml"))
json.dump({"cn_params": y["cn_params"], "cxsv_params": y["cxsv_params"], "junction_params": y["junction_params"]},
          open(out, "w"), sort_keys=True)
print(y["junction_params"]["collapse_distance"], y["junction_params"]["ins_len_ratio"])
