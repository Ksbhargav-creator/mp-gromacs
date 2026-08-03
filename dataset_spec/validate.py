#!/usr/bin/env python3
"""
validate.py -- check kernel-characterization records against the dataset schema.

Usage:
    python3 validate.py record1.json [record2.json ...]
    python3 validate.py examples/*.json

Exit code 0 if all records validate, 1 otherwise. Requires `jsonschema`
(pip install jsonschema).
"""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCHEMA_PATH = os.path.join(HERE, "kernel_characterization.schema.json")


def main(argv):
    try:
        import jsonschema
    except ImportError:
        sys.exit("jsonschema is required: pip install jsonschema")

    schema = json.load(open(SCHEMA_PATH))
    jsonschema.Draft7Validator.check_schema(schema)
    validator = jsonschema.Draft7Validator(schema)

    paths = []
    for a in argv:
        paths += sorted(glob.glob(a)) if any(c in a for c in "*?[]") else [a]
    if not paths:
        sys.exit("no record files given")

    ok = True
    for p in paths:
        try:
            rec = json.load(open(p))
        except Exception as e:
            print(f"FAIL  {p}  (not valid JSON: {e})"); ok = False; continue
        errors = sorted(validator.iter_errors(rec), key=lambda e: e.path)
        if not errors:
            print(f"OK    {p}  [{rec.get('id', '?')}]")
        else:
            ok = False
            print(f"FAIL  {p}")
            for e in errors:
                loc = "/".join(str(x) for x in e.path) or "(root)"
                print(f"        {loc}: {e.message}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
