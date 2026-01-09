# Unified Peptide Calculator

A streamlined, unified interface for calculating comprehensive physicochemical properties of peptides. This tool acts as a wrapper around the **Peptide-Tools** suite, eliminating the need for complex environment configurations or manual `PYTHONPATH` management.

## Key Features

* **Zero-Configuration:** Run immediately without setting multiple environment variables.
* **Comprehensive Property Calculation:**
* **Molecular Descriptors:** Molecular Weight (MW), LogP, Sequence Length.
* **Extinction Coefficients:** Calculated at 205 nm, 214 nm, and 280 nm.
* **Isoelectric Point (pI):** Integrated **pIChemiSt** functionality for pI and charge prediction.
* **Chemical Stability:** Automated detection of structural liabilities and formulation risks.
* **Versatile Input Support:** Accepts SMILES strings, FASTA sequences, and batch files (`.smi`, `.fasta`, `.sdf`).



## Prerequisites

Ensure the following dependencies are installed:

* **Python 3.x**
* **RDKit:** `conda install -c conda-forge rdkit` or `pip install rdkit`
* **Biopython:** `pip install biopython`

## Quick Start

### Option 1: Direct Usage (Recommended)

No installation is required. Simply navigate to the repository and run the script:

```bash
cd peptide-tools
python calculate.py --help

```

### Option 2: Installation as a Package

If you prefer to use the tool globally:

```bash
pip install -e .

```

## Usage Examples

### 1. Single Molecule Analysis

**From SMILES:**

```bash
python calculate.py --smiles "CC(C)C[C@H](NC(=O)[C@H](N)CC(O)=O)C(O)=O" --name "Asp-Leu"

```

**From FASTA:**

```bash
python calculate.py --fasta "ACDEL" --name "Peptide1"

```

### 2. Batch Processing

Process multiple sequences or structures using file inputs:

```bash
# Process a SMILES file
python calculate.py --file peptides.smi

# Process a FASTA file
python calculate.py --file sequences.fasta

# Process an SDF file
python calculate.py --file molecules.sdf

```

### 3. Output Formatting

Control how results are displayed or saved:

```bash
# Summary format (Default, human-readable)
python calculate.py --smiles "..." --format summary

# JSON format (Machine-readable)
python calculate.py --smiles "..." --format json

# Save output to a file
python calculate.py --smiles "..." --format json --output results.json

```

### 4. Advanced Configuration

Configure termini ionization for FASTA inputs:

```bash
# Treat N/C termini as capped (non-ionizable)
python calculate.py --fasta "ACDEL" --ionizable-nterm false --ionizable-cterm false

# Treat N/C termini as free (ionizable) - Default behavior
python calculate.py --fasta "ACDEL"

```

## Python API Integration

You can integrate the calculator directly into your Python scripts:

```python
from calculate import PeptideCalculator

# Initialize calculator (configure termini settings if needed)
calc = PeptideCalculator(ionizable_nterm=True, ionizable_cterm=True)

# Calculate properties from SMILES
results_smi = calc.calculate_all_properties(
    smiles="CC(C)C[C@H](NC(=O)[C@H](N)CC(O)=O)C(O)=O",
    mol_name="Asp-Leu"
)

# Calculate properties from FASTA
results_fasta = calc.calculate_all_properties(
    fasta="ACDEL",
    mol_name="Peptide1"
)

# Batch process a file
batch_results = calc.calculate_from_file("peptides.smi")

print(results_smi)

```

## File Format Specifications

* **SMILES (.smi):** Space-delimited format: `SMILES_STRING [Name]`
```text
CCNCC mol1
CC(C)O mol2

```


* **FASTA (.fasta):** Standard header and sequence format.
```text
>peptide1
ACDEL
>peptide2
FGHIK

```


* **SDF (.sdf):** Standard V2000/V3000 chemical structure file.

## Sample Output

**Summary Mode:**

```text
Molecule: Asp-Leu
  Molecular Weight: 262.263
  LogP: -2.156
  Sequence Length: 2
  FASTA: DL
  Extinction Coefficient (280nm): 0
  Extinction Coefficient (205nm): 16500
  Isoelectric Point (pI): 2.98
  Chemical Liabilities: 1 detected
    - Aspartic acid deamidation: 1 occurrence(s)

```

**JSON Mode:**

```json
{
  "mol_name": "Asp-Leu",
  "input_smiles": "CC(C)C[C@H](NC(=O)[C@H](N)CC(O)=O)C(O)=O",
  "molecular_descriptors": {
    "molecular_weight": 262.263,
    "logp": -2.156,
    "seq_length": 2,
    "fasta": "DL"
  },
  "extinction_coefficient": {
    "extn_coeff_280": 0,
    "extn_coeff_205": 16500
  },
  "pI_ChemiSt": {
    "pI": 2.98
  },
  "liabilities": {
    "Aspartic acid deamidation": {
      "liability": "Aspartic acid deamidation",
      "match_count": 1
    }
  }
}

```