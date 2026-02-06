#!/usr/bin/env python3
"""
Very small parser to extract total energy and forces from a Quantum ESPRESSO output file.
Not robust — intended as a starting point.
"""
import re
import sys

def parse_qe_energy(path):
    energy_re = re.compile(r"!    total energy\s+=\s+([-\d\.]+)")
    energies = []
    with open(path) as f:
        for line in f:
            m = energy_re.search(line)
            if m:
                energies.append(float(m.group(1)))
    return energies

def main():
    if len(sys.argv) < 2:
        print("Usage: simple_parser.py qe_output.txt")
        return
    energies = parse_qe_energy(sys.argv[1])
    for i, e in enumerate(energies):
        print(f"Step {i}: total energy = {e} Ry")

if __name__ == '__main__':
    main()