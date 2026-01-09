#!/usr/bin/env python3
"""
Unified Peptide Tools Calculator

A simplified interface to calculate all peptide properties from SMILES, FASTA, or file inputs.
No complex PYTHONPATH configuration needed.
"""

import argparse
import json
import os
import sys
from typing import Dict, Any, Optional

# Add subdirectories to path
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_SCRIPT_DIR, "pIChemiSt"))
sys.path.insert(0, os.path.join(_SCRIPT_DIR, "smi2scrambledfasta"))
sys.path.insert(0, os.path.join(_SCRIPT_DIR, "extn_coeff_fasta"))
sys.path.insert(0, os.path.join(_SCRIPT_DIR, "liabilities"))
sys.path.insert(0, os.path.join(_SCRIPT_DIR, "molecular_descriptors"))

from rdkit import Chem
from Bio import SeqIO

# Import calculation functions
from liabilities import calculate_liabilities_from_mol, calculate_liabilities_from_fasta
from extn_coeff_fasta import calc_extn_coeff_single_sequence
from molecular_descriptors import calc_molecular_descriptors_from_dict
from smi2scrambledfasta import get_scrambled_fasta_from_smiles, get_scrambled_fasta_from_mol
from pichemist.api import pichemist_from_dict
from pichemist.model import InputAttribute


class PeptideCalculator:
    """Main calculator class for all peptide properties."""
    
    def __init__(self, ionizable_nterm=True, ionizable_cterm=True):
        """
        Initialize calculator with terminal ionization settings.
        
        Args:
            ionizable_nterm: Whether N-terminus is ionizable (True) or capped (False)
            ionizable_cterm: Whether C-terminus is ionizable (True) or capped (False)
        """
        self.ionizable_nterm = ionizable_nterm
        self.ionizable_cterm = ionizable_cterm
    
    def calculate_all_properties(self, 
                                 smiles: Optional[str] = None,
                                 fasta: Optional[str] = None,
                                 mol_name: str = "Unknown") -> Dict[str, Any]:
        """
        Calculate all properties for a single peptide.
        
        Args:
            smiles: SMILES string (optional)
            fasta: FASTA sequence (optional)
            mol_name: Name of the molecule
            
        Returns:
            Dictionary containing all calculated properties
        """
        if not smiles and not fasta:
            raise ValueError("Either SMILES or FASTA must be provided")
        
        results = {
            "mol_name": mol_name,
            "input_smiles": smiles,
            "input_fasta": fasta,
        }
        
        # Create mol_supply_json format for compatibility
        mol_supply_json = self._create_mol_supply_json(smiles, fasta, mol_name)
        
        # Calculate molecular descriptors
        try:
            descriptors = calc_molecular_descriptors_from_dict(mol_supply_json)
            results["molecular_descriptors"] = descriptors.get(1, {})
        except Exception as e:
            results["molecular_descriptors"] = {"error": str(e)}
        
        # Calculate extinction coefficient (requires FASTA)
        try:
            if fasta:
                ec_dict = calc_extn_coeff_single_sequence(fasta)
            elif smiles:
                # Convert SMILES to FASTA first
                derived_fasta = results["molecular_descriptors"].get("fasta", "")
                if derived_fasta:
                    ec_dict = calc_extn_coeff_single_sequence(derived_fasta)
                else:
                    ec_dict = {"error": "Could not derive FASTA from SMILES"}
            results["extinction_coefficient"] = ec_dict
        except Exception as e:
            results["extinction_coefficient"] = {"error": str(e)}
        
        # Calculate pI (isoelectric point)
        try:
            pi_dict = pichemist_from_dict(
                mol_supply_json,
                method="pkamatcher",
                ionizable_nterm=self.ionizable_nterm,
                ionizable_cterm=self.ionizable_cterm,
                print_fragments=False,
                plot_ph_q_curve=False,
                generate_fragment_images=False
            )
            results["pI_ChemiSt"] = pi_dict.get(1, {})
        except Exception as e:
            results["pI_ChemiSt"] = {"error": str(e)}
        
        # Calculate liabilities
        try:
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                liabilities = calculate_liabilities_from_mol(mol)
            elif fasta:
                liabilities = calculate_liabilities_from_fasta(fasta)
            results["liabilities"] = liabilities
        except Exception as e:
            results["liabilities"] = {"error": str(e)}
        
        return results
    
    def _create_mol_supply_json(self, smiles, fasta, mol_name):
        """Create the mol_supply_json format expected by internal functions."""
        mol_obj = None
        derived_fasta = fasta
        
        if smiles:
            mol_obj = Chem.MolFromSmiles(smiles)
            if mol_obj is None:
                raise ValueError(f"Invalid SMILES: {smiles}")
            if not fasta:
                derived_fasta = get_scrambled_fasta_from_mol(mol_obj)
        
        return {
            1: {
                "mol_name": mol_name,
                "mol_obj": mol_obj,
                "fasta": derived_fasta
            }
        }
    
    def calculate_from_file(self, filepath: str) -> Dict[int, Dict[str, Any]]:
        """
        Calculate properties for all molecules in a file.
        
        Args:
            filepath: Path to input file (.smi, .fasta, or .sdf)
            
        Returns:
            Dictionary mapping molecule index to properties
        """
        ext = os.path.splitext(filepath)[1].lower()
        
        if ext == ".smi":
            return self._process_smi_file(filepath)
        elif ext == ".fasta":
            return self._process_fasta_file(filepath)
        elif ext == ".sdf":
            return self._process_sdf_file(filepath)
        else:
            raise ValueError(f"Unsupported file format: {ext}")
    
    def _process_smi_file(self, filepath: str) -> Dict[int, Dict[str, Any]]:
        """Process SMILES file."""
        results = {}
        with open(filepath, 'r') as f:
            for idx, line in enumerate(f, start=1):
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split()
                smiles = parts[0]
                mol_name = parts[1] if len(parts) > 1 else f"mol_{idx}"
                
                try:
                    results[idx] = self.calculate_all_properties(
                        smiles=smiles, 
                        mol_name=mol_name
                    )
                except Exception as e:
                    results[idx] = {
                        "mol_name": mol_name,
                        "error": str(e)
                    }
        return results
    
    def _process_fasta_file(self, filepath: str) -> Dict[int, Dict[str, Any]]:
        """Process FASTA file."""
        results = {}
        with open(filepath, 'r') as f:
            for idx, record in enumerate(SeqIO.parse(f, "fasta"), start=1):
                fasta_seq = str(record.seq)
                mol_name = record.id
                
                try:
                    results[idx] = self.calculate_all_properties(
                        fasta=fasta_seq,
                        mol_name=mol_name
                    )
                except Exception as e:
                    results[idx] = {
                        "mol_name": mol_name,
                        "error": str(e)
                    }
        return results
    
    def _process_sdf_file(self, filepath: str) -> Dict[int, Dict[str, Any]]:
        """Process SDF file."""
        results = {}
        supplier = Chem.SDMolSupplier(filepath)
        
        for idx, mol in enumerate(supplier, start=1):
            if mol is None:
                results[idx] = {
                    "mol_name": f"mol_{idx}",
                    "error": "Could not parse molecule"
                }
                continue
            
            mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{idx}"
            smiles = Chem.MolToSmiles(mol)
            
            try:
                results[idx] = self.calculate_all_properties(
                    smiles=smiles,
                    mol_name=mol_name
                )
            except Exception as e:
                results[idx] = {
                    "mol_name": mol_name,
                    "error": str(e)
                }
        
        return results


def print_results(results: Dict[str, Any], format: str = "json"):
    """Print results in specified format."""
    if format == "json":
        print(json.dumps(results, indent=2, default=str))
    elif format == "summary":
        if "mol_name" in results:
            # Single molecule
            _print_single_summary(results)
        else:
            # Multiple molecules
            for idx, mol_results in results.items():
                print(f"\n{'='*80}")
                print(f"Molecule {idx}")
                print(f"{'='*80}")
                _print_single_summary(mol_results)


def _print_single_summary(results: Dict[str, Any]):
    """Print summary for a single molecule."""
    print(f"\nMolecule: {results.get('mol_name', 'Unknown')}")
    
    if "error" in results:
        print(f"  ERROR: {results['error']}")
        return
    
    # Molecular descriptors
    desc = results.get("molecular_descriptors", {})
    if "error" not in desc:
        print(f"  Molecular Weight: {desc.get('molecular_weight', 'N/A')}")
        print(f"  LogP: {desc.get('logp', 'N/A')}")
        print(f"  Sequence Length: {desc.get('seq_length', 'N/A')}")
        print(f"  FASTA: {desc.get('fasta', 'N/A')}")
    
    # Extinction coefficient
    ec = results.get("extinction_coefficient", {})
    if "error" not in ec:
        print(f"  Extinction Coefficient (280nm): {ec.get('extn_coeff_280', 'N/A')}")
        print(f"  Extinction Coefficient (205nm): {ec.get('extn_coeff_205', 'N/A')}")
    
    # pI
    pi = results.get("pI_ChemiSt", {})
    if "error" not in pi:
        print(f"  Isoelectric Point (pI): {pi.get('pI', 'N/A')}")
    
    # Liabilities
    liabilities = results.get("liabilities", {})
    if liabilities and "error" not in liabilities:
        print(f"  Chemical Liabilities: {len(liabilities)} detected")
        for liability_name, liability_info in liabilities.items():
            print(f"    - {liability_name}: {liability_info.get('match_count', 0)} occurrence(s)")


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Calculate all physicochemical properties for peptides",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Calculate from SMILES
  python calculate.py --smiles "CC(C)C[C@H](NC(=O)[C@H](N)CC(O)=O)C(O)=O" --name "Asp-Leu"
  
  # Calculate from FASTA
  python calculate.py --fasta "ACDEL" --name "Peptide1"
  
  # Calculate from file
  python calculate.py --file peptides.smi
  python calculate.py --file sequences.fasta
  
  # Output format options
  python calculate.py --smiles "..." --format json
  python calculate.py --smiles "..." --format summary
        """
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--smiles", type=str, help="Input SMILES string")
    input_group.add_argument("--fasta", type=str, help="Input FASTA sequence")
    input_group.add_argument("--file", type=str, help="Input file (.smi, .fasta, or .sdf)")
    
    parser.add_argument("--name", type=str, default="Unknown", 
                       help="Molecule name (for single inputs)")
    parser.add_argument("--format", type=str, choices=["json", "summary"], 
                       default="summary", help="Output format")
    parser.add_argument("--output", type=str, help="Output file (optional, prints to stdout if not specified)")
    parser.add_argument("--ionizable-nterm", type=lambda x: x.lower() == "true",
                       default=True, help="N-terminus is ionizable (true/false)")
    parser.add_argument("--ionizable-cterm", type=lambda x: x.lower() == "true",
                       default=True, help="C-terminus is ionizable (true/false)")
    
    args = parser.parse_args()
    
    # Initialize calculator
    calculator = PeptideCalculator(
        ionizable_nterm=args.ionizable_nterm,
        ionizable_cterm=args.ionizable_cterm
    )
    
    # Calculate properties
    try:
        if args.file:
            results = calculator.calculate_from_file(args.file)
        elif args.smiles:
            results = calculator.calculate_all_properties(
                smiles=args.smiles,
                mol_name=args.name
            )
        elif args.fasta:
            results = calculator.calculate_all_properties(
                fasta=args.fasta,
                mol_name=args.name
            )
        
        # Output results
        if args.output:
            with open(args.output, 'w') as f:
                if args.format == "json":
                    json.dump(results, f, indent=2, default=str)
                else:
                    # Redirect stdout to file temporarily
                    old_stdout = sys.stdout
                    sys.stdout = f
                    print_results(results, format=args.format)
                    sys.stdout = old_stdout
            print(f"Results written to {args.output}")
        else:
            print_results(results, format=args.format)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
