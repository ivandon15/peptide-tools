#!/usr/bin/env python3
"""
Quick start examples for the unified peptide calculator.
Run these examples to see how to use the calculator.
"""

from calculate import PeptideCalculator
import json


def example_1_simple_smiles():
    """Example 1: Calculate properties from SMILES"""
    print("=" * 80)
    print("Example 1: Calculate from SMILES")
    print("=" * 80)
    
    calc = PeptideCalculator()
    
    # Dipeptide: Asp-Leu
    smiles = "CC(C)C[C@H](NC(=O)[C@H](N)CC(O)=O)C(O)=O"
    
    results = calc.calculate_all_properties(
        smiles=smiles,
        mol_name="Asp-Leu"
    )
    
    print("\nResults:")
    print(json.dumps(results, indent=2, default=str))


def example_2_simple_fasta():
    """Example 2: Calculate properties from FASTA"""
    print("\n" + "=" * 80)
    print("Example 2: Calculate from FASTA")
    print("=" * 80)
    
    calc = PeptideCalculator()
    
    results = calc.calculate_all_properties(
        fasta="ACDEL",
        mol_name="Pentapeptide"
    )
    
    print("\nResults:")
    print(json.dumps(results, indent=2, default=str))


def example_3_batch_processing():
    """Example 3: Batch process multiple peptides"""
    print("\n" + "=" * 80)
    print("Example 3: Batch Processing")
    print("=" * 80)
    
    calc = PeptideCalculator()
    
    peptides = [
        ("ACDEL", "Peptide1"),
        ("FGHIK", "Peptide2"),
        ("LMNPR", "Peptide3"),
    ]
    
    all_results = {}
    for idx, (fasta, name) in enumerate(peptides, start=1):
        results = calc.calculate_all_properties(fasta=fasta, mol_name=name)
        all_results[idx] = results
        
        # Print summary
        print(f"\n{name} ({fasta}):")
        print(f"  MW: {results['molecular_descriptors']['molecular_weight']}")
        print(f"  pI: {results['pI_ChemiSt']['pI']}")
        print(f"  Liabilities: {len(results['liabilities'])}")


def example_4_extract_specific_properties():
    """Example 4: Extract specific properties you need"""
    print("\n" + "=" * 80)
    print("Example 4: Extract Specific Properties")
    print("=" * 80)
    
    calc = PeptideCalculator()
    
    fasta = "PEPTIDE"
    results = calc.calculate_all_properties(fasta=fasta, mol_name="MyPeptide")
    
    # Extract only what you need
    mol_weight = results['molecular_descriptors']['molecular_weight']
    pi = results['pI_ChemiSt']['pI']
    ec_280 = results['extinction_coefficient']['extn_coeff_280']
    logp = results['molecular_descriptors']['logp']
    
    print(f"\nExtracted properties for {fasta}:")
    print(f"  Molecular Weight: {mol_weight} Da")
    print(f"  Isoelectric Point: {pi}")
    print(f"  Extinction Coefficient (280nm): {ec_280} M⁻¹cm⁻¹")
    print(f"  LogP: {logp}")


def example_5_terminal_capping():
    """Example 5: Compare free vs capped terminals"""
    print("\n" + "=" * 80)
    print("Example 5: Terminal Capping Effects")
    print("=" * 80)
    
    fasta = "PEPTIDE"
    
    # Free terminals (default)
    calc_free = PeptideCalculator(ionizable_nterm=True, ionizable_cterm=True)
    results_free = calc_free.calculate_all_properties(fasta=fasta, mol_name="Free")
    
    # Capped terminals
    calc_capped = PeptideCalculator(ionizable_nterm=False, ionizable_cterm=False)
    results_capped = calc_capped.calculate_all_properties(fasta=fasta, mol_name="Capped")
    
    print(f"\nComparison for {fasta}:")
    print(f"  Free terminals:")
    print(f"    pI: {results_free['pI_ChemiSt']['pI']}")
    print(f"    MW: {results_free['molecular_descriptors']['molecular_weight']}")
    
    print(f"  Capped terminals:")
    print(f"    pI: {results_capped['pI_ChemiSt']['pI']}")
    print(f"    MW: {results_capped['molecular_descriptors']['molecular_weight']}")


def example_6_liability_details():
    """Example 6: Analyze chemical liabilities in detail"""
    print("\n" + "=" * 80)
    print("Example 6: Chemical Liability Analysis")
    print("=" * 80)
    
    calc = PeptideCalculator()
    
    # Sequence with known liabilities
    fasta = "DGNGM"  # Contains D-G (Asp deamidation), N-G (Asn deamidation), M (oxidation)
    results = calc.calculate_all_properties(fasta=fasta, mol_name="Unstable_Peptide")
    
    print(f"\nLiability analysis for {fasta}:")
    liabilities = results['liabilities']
    
    if liabilities:
        print(f"  Found {len(liabilities)} potential liabilities:")
        for liability_name, liability_info in liabilities.items():
            print(f"\n  {liability_name}:")
            print(f"    Pattern: {liability_info['fasta_pattern']}")
            print(f"    Count: {liability_info['match_count']}")
            print(f"    Stability factors: {liability_info['factors_favoring_stability']}")
            print(f"    Instability factors: {liability_info['factors_disfavoring_stability']}")
    else:
        print("  No liabilities detected")


if __name__ == "__main__":
    print("\n" + "=" * 80)
    print("UNIFIED PEPTIDE CALCULATOR - QUICK START EXAMPLES")
    print("=" * 80 + "\n")
    
    examples = [
        example_1_simple_smiles,
        example_2_simple_fasta,
        example_3_batch_processing,
        example_4_extract_specific_properties,
        example_5_terminal_capping,
        example_6_liability_details,
    ]
    
    for example in examples:
        try:
            example()
        except Exception as e:
            print(f"\n✗ Example failed: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 80)
    print("Examples completed!")
    print("=" * 80)
