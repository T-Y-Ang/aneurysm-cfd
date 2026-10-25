"""
Correct the OpenFOAM vessel-wall boundary classification.

Gmsh mesh imports may classify vessel_wall as a generic 'patch'.
OpenFOAM's wallShearStress function requires a wall-type mesh patch.

This utility changes the classification without altering mesh geometry,
face counts, face ordering, or other boundary regions.
"""

from pathlib import Path
import re
import shutil
import sys


# Default location of the optimised patient-derived CFD case.
BOUNDARY_FILE = Path(
    "cases/patient_aneurysm_optimized/constant/polyMesh/boundary"
)


def main():
    # Confirm that the expected boundary file exists.
    if not BOUNDARY_FILE.is_file():
        sys.exit(f"ERROR: Boundary file not found: {BOUNDARY_FILE}")

    original = BOUNDARY_FILE.read_text()

    # Find only the vessel_wall dictionary.
    pattern = re.compile(
        r"(\bvessel_wall\s*\{)([^{}]*)(\})",
        re.MULTILINE,
    )

    matches = list(pattern.finditer(original))

    if len(matches) != 1:
        sys.exit(
            f"ERROR: Expected one vessel_wall definition; found {len(matches)}."
        )

    match = matches[0]
    block = match.group(2)

    # Require the expected patch classification before modifying it.
    type_match = re.search(r"(?m)^(\s*type\s+)(\w+)(\s*;)", block)

    if type_match is None:
        sys.exit("ERROR: vessel_wall has no recognised type entry.")

    current_type = type_match.group(2)

    if current_type == "wall":
        print("vessel_wall is already classified as wall. No change needed.")
        return

    if current_type != "patch":
        sys.exit(f"ERROR: Unexpected vessel_wall type: {current_type}")

    # Replace only the patch type within vessel_wall.
    updated_block = re.sub(
        r"(?m)^(\s*type\s+)patch(\s*;)",
        r"\1wall\2",
        block,
        count=1,
    )

    # Update physicalType if present.
    updated_block = re.sub(
        r"(?m)^(\s*physicalType\s+)patch(\s*;)",
        r"\1wall\2",
        updated_block,
        count=1,
    )

    # Preserve everything outside the vessel_wall dictionary.
    updated = (
        original[:match.start(2)]
        + updated_block
        + original[match.end(2):]
    )

    # Preserve the original file before writing.
    backup = BOUNDARY_FILE.with_name("boundary.backup")
    shutil.copy2(BOUNDARY_FILE, backup)

    BOUNDARY_FILE.write_text(updated)

    print("Updated vessel_wall classification: patch -> wall")
    print(f"Backup saved to: {backup}")


if __name__ == "__main__":
    main()