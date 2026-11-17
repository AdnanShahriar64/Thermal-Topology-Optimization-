#!/usr/bin/env python3
"""Insert pre-solver ground radiation preparation into an existing Runner_4.m.

Leaves the original MATLAB solver (starting at %% Solver) and all existing
functions, including Surface_Integration_3D_Flux_Thermal, byte-for-byte
unchanged. Place the supplied Ground_Radiation_*.m helpers on MATLAB's path.

Usage:
    python apply_ground_radiation_patch.py /path/to/Runner_4.m
    python apply_ground_radiation_patch.py /path/to/Thermal_Model_Runner_Complete.m

Writes a new file beside the original with '_ground_radiation' in its name.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import re

ANCHOR = '    [K_Conv] = Assembly(Ke_Convection,1,Mesh_Information,Element_Information);'
BEGIN = '%% Ground IR and ground-reflected solar (ALBEDO) - PRE-SOLVER ONLY'


def patch(source: Path, destination: Path, block: Path) -> None:
    raw = source.read_bytes()
    try:
        contents = raw.decode('utf-8')
    except UnicodeDecodeError:
        contents = raw.decode('cp1252')
    if BEGIN in contents:
        raise ValueError('The ground-radiation setup already exists. Refusing to patch twice.')
    if contents.count(ANCHOR) != 1:
        raise ValueError('Expected exactly one K_Conv assembly anchor; found %d.' % contents.count(ANCHOR))
    solver_re = re.compile(r'(?m)^[ \t]*%% Solver[ \t]*$')
    original_solver = solver_re.search(contents, contents.index(ANCHOR))
    if original_solver is None:
        raise ValueError('Could not find a standalone %% Solver heading after the assembly block.')
    original_solver_and_beyond = contents[original_solver.start():]
    insertion = block.read_text(encoding='utf-8').strip('\n')
    newline = '\r\n' if b'\r\n' in raw else '\n'
    insertion = insertion.replace('\n', newline)
    modified = contents.replace(ANCHOR, ANCHOR + newline + newline + insertion, 1)
    modified_solver = solver_re.search(modified, modified.index(ANCHOR))
    assert modified_solver is not None
    assert modified[modified_solver.start():] == original_solver_and_beyond
    if destination.resolve() == source.resolve():
        raise ValueError('Destination must differ from input; the original is preserved.')
    destination.write_bytes(modified.encode('utf-8'))
    print(f'Created {destination}')
    print('Verified: original %% Solver and everything after it is unchanged.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path, help='Path to existing MATLAB .m file')
    parser.add_argument('--output', type=Path, help='Output .m file (default: *_ground_radiation.m)')
    args = parser.parse_args()
    source = args.source
    if not source.is_file():
        parser.error(f'Source file not found: {source}')
    destination = args.output or source.with_name(source.stem + '_ground_radiation.m')
    patch(source, destination, Path(__file__).with_name('Thermal_Model_Runner_Complete_presolver_insert.txt'))


if __name__ == '__main__':
    main()
