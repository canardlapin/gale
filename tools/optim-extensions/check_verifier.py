#!/usr/bin/env python3
"""Reject corrupted certificates in retained fixture results (no third-party packages)."""
import csv
import tempfile
from pathlib import Path
from compare import OUT, verify


def main():
    with (OUT / "jvm.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    with tempfile.TemporaryDirectory(prefix="gale-optim-verifier-") as folder:
        for field, value in [("residual", "NaN"), ("residual", "1"), ("status", "Unknown"), ("objective", "NaN")]:
            altered = [dict(row) for row in rows]
            altered[0][field] = value
            path = Path(folder) / "corrupt.tsv"
            with path.open("w") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(altered)
            try:
                verify(path, OUT / "python.json", Path(folder) / "result")
            except ValueError:
                print(f"Rejected {field}={value}")
            else:
                raise AssertionError(f"accepted corrupt {field}={value}")
    print("All four verifier corruption probes pass")


if __name__ == "__main__":
    main()
