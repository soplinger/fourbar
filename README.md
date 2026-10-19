# Fourbar

**Offline four-bar linkage analyser: describe the links in a plain-text file and get the mechanism class, crank and rocker ranges, transmission angles, dead points, coupler-point paths, and SVG drawings.**

Fourbar reads a small `.lnk` file, sweeps the crank through a full turn, and tells you what the linkage really does: whether the crank turns fully, how far the rocker swings, how close it gets to a toggle that locks it, and where any coupler point travels. It is a single Python file using only the standard library. No servers, accounts, API keys, or network calls.

## Features

- **Plain-text input** with `#` comments, a `name`, optional `units`, and four link lengths.
- **Line-numbered lint**: every problem in a file is reported at once as `file:line: message`, including unknown directives, duplicates, bad numbers, reserved point names, and missing links.
- **Grashof test and mechanism class** (`check`, `analyse`): crank-rocker, rocker-crank, double-crank, double-rocker, or triple-rocker, with a change-point flag.
- **Joint positions at any crank angle** (`at`): O2, A, B, O4, rocker angle, transmission angle, and every coupler point.
- **Assembly branch choice** (`--branch open|crossed`) on every analysis command.
- **Crank ranges of motion** from a full-turn feasibility sweep, reported as angle ranges.
- **Rocker swing or full rotation** (`analyse`), computed on an unwrapped angle so there are no 360° jumps.
- **Transmission angle extremes** with the angle where each occurs, checked against the file's `min_transmission` limit.
- **Near-toggle (dead point) detection** at the crank angles where the transmission angle nearly reaches 0° or 180°.
- **Coupler points** defined by distance and angle along the coupler, with bounding boxes over the sweep (`analyse`).
- **ASCII path plot** of any coupler point (`trace`).
- **SVG drawing** of the mechanism at one angle, optionally with the traced path of a point (`svg`, `-o FILE`).
- **CSV sweep export** of every sampled angle for spreadsheets or plotting (`sweep`).
- **JSON analysis** for scripts (`analyse --json`).
- **CI gate** that fails the build when the minimum transmission angle is below a limit (`gate`, `--min-transmission`), with documented exit codes.

## Quick start

Write a linkage file. The bundled example is `examples/crank-rocker.lnk`:

```text
name   wiper-linkage
units  mm

ground   4.0    # fixed distance between O2 and O4
crank    1.0    # O2 to A, the input link
coupler  3.0    # A to B
rocker   2.5    # O4 to B, the output link

min_transmission 40

point  T  3.0   0     # coupler tip, sits on B
point  P  1.5   30    # a point off the coupler line
```

Check it and see the mechanism class:

```sh
python3 fourbar.py check examples/crank-rocker.lnk
```

```text
OK examples/crank-rocker.lnk: wiper-linkage is a crank-rocker (crank rotates fully, rocker oscillates)
```

Run the full analysis:

```sh
python3 fourbar.py analyse examples/crank-rocker.lnk
```

It prints the Grashof sums, the crank range, the rocker swing, the transmission-angle extremes against the 40° limit, any near-toggle positions, and the bounding box of each coupler point.

Other useful commands:

```sh
python3 fourbar.py at examples/crank-rocker.lnk 90          # joints at crank = 90 deg
python3 fourbar.py trace examples/crank-rocker.lnk P        # ASCII path of point P
python3 fourbar.py svg examples/crank-rocker.lnk 90 --trace P -o linkage.svg
python3 fourbar.py sweep examples/crank-rocker.lnk > sweep.csv
python3 fourbar.py gate examples/crank-rocker.lnk --min-transmission 45
```

## Linkage file format

| Directive | Meaning |
| --- | --- |
| `name TEXT` | Display name (optional). |
| `units WORD` | Length unit label, for display only (optional). |
| `ground N` | Distance between the fixed pivots O2 and O4. Required. |
| `crank N` | Input link O2A. Required. |
| `coupler N` | Floating link AB. Required. |
| `rocker N` | Output link O4B. Required. |
| `min_transmission DEG` | Gate limit in degrees, in (0, 90]. Default 40. |
| `point NAME DIST ANGLE` | Coupler point at DIST along AB from A, rotated ANGLE degrees from AB. Names `A`, `B`, `O2`, `O4` are reserved. |

All lengths must be positive finite numbers. Use one unit system throughout.

## Command reference

| Command | Purpose | Exit codes |
| --- | --- | --- |
| `check FILE` | Lint and print the mechanism class | 0 ok, 1 file error |
| `analyse FILE [--json]` | Full report, optionally as JSON | 0 ok, 1 file error |
| `at FILE ANGLE` | Joint positions at one crank angle | 0 ok, 1 file error, 2 links cannot close |
| `sweep FILE` | CSV of every sampled crank angle | 0 ok, 1 file error |
| `trace FILE POINT` | ASCII plot of one coupler point | 0 ok, 1 file or point error, 2 nothing to plot |
| `svg FILE ANGLE [--trace POINT] [-o OUT]` | SVG drawing | 0 ok, 1 file or point error, 2 links cannot close |
| `gate FILE [--min-transmission DEG]` | CI gate on the minimum transmission angle | 0 pass, 1 file error, 2 fail |

Shared options for `analyse`, `at`, `sweep`, `trace`, `svg`, and `gate`:

- `--step DEG` sets the sweep resolution, between 0.05 and 90 (default 1).
- `--branch open|crossed` selects the assembly branch (default open).

Usage errors such as a bad `--step` exit with code 2 from argparse.

## How it works

- The fixed pivots are O2 = (0, 0) and O4 = (ground, 0). The crank angle θ places A.
- The coupler circle about A and the rocker circle about O4 intersect at B. Solving gives the rocker angle φ. The `open` branch takes the `+` solution and `crossed` the `-` solution.
- The transmission angle μ is the angle at B between BA and BO4. Its cosine follows from the law of cosines, so it depends only on the distance AO4.
- A crank angle is feasible when the two circles intersect. The sweep samples evenly over 360°, so 0°, 90°, 180°, and 270° are sampled exactly.
- Rocker swing is measured on an unwrapped angle for each feasible run of crank angles.

## Limitations

- The input is always the crank. For a rocker-crank mechanism the crank oscillates, so the feasible crank range is partial and the rocker's full rotation is reported by the class, not by the sweep.
- Positions are sampled, not solved exactly, so dead points are reported as near-toggle positions within 1° and the sweep resolution.
- Only the four-bar linkage is modelled. Coupler points and drawings are for visual and kinematic checks, not for dynamic force analysis.

## Testing

The test suite uses pytest and checks parsing, lint messages with line numbers, Grashof classes, joint geometry at many angles, the transmission angle against the law of cosines, feasibility limits, dead-point detection, the CLI exit codes, and the output formats.

```sh
python3 -m pytest -q
```

Run it from the project root, or from anywhere since the tests locate the module themselves.

## Project layout

```text
fourbar.py                  the analyser and CLI (single file)
examples/crank-rocker.lnk   a valid demo linkage
examples/broken.lnk         a deliberately broken file used by the tests
tests/test_fourbar.py       test suite
README.md                   this file
```
