#!/usr/bin/env python3
"""
Simple SCF convergence plotter.
Reads a plain text file with two columns: iteration energy_change
and plots energy vs iteration and residual vs iteration.

Usage:
  python scf_convergence_plotter.py scf_log.txt
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_data(path):
    data = np.loadtxt(path)
    iterations = np.arange(1, data.shape[0]+1)
    energy = data[:,0]
    resid = data[:,1] if data.shape[1]>1 else None
    return iterations, energy, resid

def main():
    if len(sys.argv) < 2:
        print("Usage: python scf_convergence_plotter.py scf_data.txt")
        sys.exit(1)
    it, en, r = read_data(sys.argv[1])
    plt.figure()
    plt.plot(it, en, marker='o')
    plt.xlabel('Iteration')
    plt.ylabel('Energy (a.u.) or ΔE')
    plt.title('SCF energy convergence')
    plt.grid(True)
    if r is not None:
        plt.twinx().plot(it, r, 'r--', label='residual')
    plt.show()

if __name__ == '__main__':
    main()