#!/usr/bin/env python3
"""
Setup script for the unified peptide calculator.
This installs all required dependencies.
"""

from setuptools import setup, find_packages

setup(
    name="peptide-calculator",
    version="1.0.0",
    description="Unified calculator for peptide physicochemical properties",
    author="Kelsey",
    py_modules=["calculate"],
    install_requires=[
        "rdkit>=2022.3.1",
        "biopython>=1.79",
    ],
    entry_points={
        "console_scripts": [
            "peptide-calc=calculate:main",
        ],
    },
    python_requires=">=3.7",
)
