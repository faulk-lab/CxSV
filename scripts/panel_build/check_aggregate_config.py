#!/usr/bin/env python3
"""Check the copy-number run-configuration keys of a rebuild's aggregate.yaml
(rebuild_downstream.sh, before any Snakemake run).

    check_aggregate_config.py AGGREGATE_YAML

aggregate.yaml is JSON (make_aggregate_yaml.py). cn_depth_precomputed (a directory)
and cn_strata_table (a file) must be set, not "none", and exist. Exit 1 with one
line per problem otherwise.
"""
import json
import os
import sys

REQUIRED = (("cn_depth_precomputed", os.path.isdir), ("cn_strata_table", os.path.isfile))


def problems(path):
    try:
        cfg = json.load(open(path))
    except (OSError, ValueError) as e:
        return [f"{path}: not readable as JSON ({e})"]
    out = []
    for key, ok in REQUIRED:
        v = cfg.get(key)
        if not v or str(v).lower() == "none":
            out.append(f"{path}: {key} unset")
        elif not ok(v):
            out.append(f"{path}: {key} = {v} does not exist")
    return out


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        sys.exit(__doc__)
    p = problems(argv[0])
    for line in p:
        print(line, file=sys.stderr)
    return 1 if p else 0


if __name__ == "__main__":
    sys.exit(main())
