#!/usr/bin/env python3
"""cn_params.json of the integration branch (as rule cn_params_json) + comparison with the release's."""
import json, sys, yaml
code, out, rel = sys.argv[1:4]
y = yaml.safe_load(open(f"{code}/config/calling.yaml"))
body = {"cn_params": y["cn_params"], "cxsv_params": y["cxsv_params"], "junction_params": y["junction_params"]}
json.dump(body, open(out, "w"), sort_keys=True)
r = json.load(open(rel))
def flat(d, p=""):
    o = {}
    for k, v in d.items():
        if isinstance(v, dict): o.update(flat(v, p + k + "."))
        else: o[p + k] = v
    return o
a, b = flat(r), flat(body)
print("only in release:", sorted(set(a) - set(b))); print("only in branch:", sorted(set(b) - set(a)))
print("differ:", {k: (a[k], b[k]) for k in set(a) & set(b) if a[k] != b[k]})
