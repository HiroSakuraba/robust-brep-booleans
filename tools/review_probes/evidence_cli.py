#!/usr/bin/env python3
"""Display a brepkernel evidence record (G17).

Shows the verdict, certification block, and input hashes.  The
conspicuous "unprobed": true marker is always surfaced when present.

Usage:
    python tools/review_probes/evidence_cli.py <record.json>
    python tools/review_probes/evidence_cli.py --validate <record.json>
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from brepkernel import evidence


def show(record):
    print(f"schema:  {record.get('schema')}")
    print(f"name:    {record.get('name')}")
    out = record.get("outcome", {})
    print(f"verdict: {out.get('verdict')}")
    if out.get("verdict") == "refused":
        print(f"  category: {out.get('category')}")
        print(f"  stage:    {out.get('stage')}")
    cert = record.get("certification", {})
    print("certification:")
    print(f"  mode:               {cert.get('mode')}")
    print(f"  completeness_probe: {cert.get('completeness_probe')}")
    print(f"  allow_nonmanifold:  {cert.get('allow_nonmanifold')}")
    if cert.get("unprobed"):
        print("  *** UNPROBED: completeness probe did not run ***")
    for blk in record.get("inputs", []):
        print(f"input {blk.get('role')}: {blk.get('sha256', '')[:16]}... "
              f"({blk.get('brep_bytes')} bytes)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("record", help="path to evidence JSON record")
    ap.add_argument("--validate", action="store_true",
                    help="run validate_evidence and report problems")
    args = ap.parse_args()
    with open(args.record) as f:
        record = json.load(f)
    show(record)
    if args.validate:
        problems = evidence.validate_evidence(record)
        if problems:
            print("\nvalidation problems:")
            for p in problems:
                print(f"  - {p}")
            sys.exit(1)
        print("\nrecord validates clean")


if __name__ == "__main__":
    main()
