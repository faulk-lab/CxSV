#!/usr/bin/env bash
# Print "<tool>\t<version>" for every tool whose output defines the calls.
# Part of the calling fingerprint (workflow/lib/cxsv_fingerprint.py): a panel
# and a patient run must report identical lines.
set -uo pipefail
v() { printf '%s\t%s\n' "$1" "$2"; }
v sniffles "$(sniffles --version 2>&1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)"
v cuteSV   "$(cuteSV --version 2>&1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)"
v svim     "$(svim --version 2>&1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)"
v SURVIVOR "$(SURVIVOR 2>&1 | grep -i 'version' | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)"
v bcftools "$(bcftools --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)"
v bedtools "$(bedtools --version 2>&1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)"
v samtools "$(samtools --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)"
