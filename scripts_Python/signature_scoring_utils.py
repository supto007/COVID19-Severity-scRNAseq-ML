import pandas as pd
import numpy as np

def calculate_covid_severity_score(expression_df):
    """
    Independent function to calculate the S100A8/S100A9/ITGB2 Signature Score.
    Uses Z-score Stacking Logic.
    """
    target_genes = ['S100A8', 'S100A9', 'ITGB2']
    
    # Check gene availability
    available_genes = [g for g in target_genes if g in expression_df.columns]
    
    if len(available_genes) < 3:
        print(f"Warning: Only {available_genes} found in dataset.")
    
    # Step 1: Subset and Z-score normalization
    subset = expression_df[available_genes]
    z_scored = (subset - subset.mean()) / subset.std()
    
    # Step 2: Composite Score calculation (Stacking)
    signature_score = z_scored.mean(axis=1)
    
    return signature_score

print("Utility functions for Signature Scoring are ready.")
