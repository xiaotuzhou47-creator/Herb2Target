# -*- coding: utf-8 -*-
"""Herb2Target Validation Script
Run: python run_validation.py
Records: 2026-06-15 - First full re-run after paper revision
"""

import sys, os, json, math, warnings
sys.stdout.reconfigure(encoding='utf-8')
os.environ['FLASK_ENV'] = 'development'
warnings.filterwarnings('ignore')

# Change to script directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Add stubs for missing modules before importing app
import types
# uniprot_pdb stub
uniprot_pdb_mod = types.ModuleType('uniprot_pdb')
uniprot_pdb_mod.UNIPROT_PDB = {}
sys.modules['uniprot_pdb'] = uniprot_pdb_mod

# Import the app's prediction function
sys.path.insert(0, '.')
from app import predict_molecule_targets

# ======== 10 Validation Pairs from Table S6 ========
validation_pairs = [
    ('Resveratrol', 'C1=CC(=CC=C1C=CC2=CC(=CC(=C2)O)O)O', 'PTGS2'),
    ('Curcumin', 'COC1=C(C=CC(=C1)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC)O', 'TNF'),
    ('Berberine', 'COC1=C(C2=C(C=C1)C=C3C4=CC5=C(C=C4CC[N+]3=C2)OCO5)OC', 'AKT1'),
    ('Baicalein', 'C1=CC=C(C=C1)C2=CC(=O)C3=C(C=C(C(=C3O2)O)O)O', 'PTGS2'),
    ('Kaempferol', 'C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O', 'AKT1'),
    ('Hesperetin', 'COC1=C(C=CC(=C1)C2CC(=O)C3=C(C=C(C=C3O2)O)O)O', 'PTGS2'),
    ('Apigenin', 'C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O', 'PTGS2'),
    ('Wogonin', 'COC1=C(C=C(C2=C1OC(=CC2=O)C3=CC=CC=C3)O)O', 'PTGS2'),
    ('EGCG', 'C1=C(C=C(C(=C1O)O)O)C2CC(=O)C3=C(C=C(C(=C3O2)O)O)O', 'EGFR'),
    ('Luteolin', 'C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O)O', 'TNF'),
]

print()
print('='*70)
print('HERB2TARGET VALIDATION RUN')
print('Date: 2026-06-15')
print('Source: /e/demo/demo/app.py predict_molecule_targets()')
print('='*70)
print()

ht_ranks = []
for name, smi, expected_target in validation_pairs:
    try:
        result = predict_molecule_targets(smi, top_n=30)
        targets = result.get('targets', [])
        target_names = [t['gene_name'] for t in targets]

        if expected_target in target_names:
            rank = target_names.index(expected_target) + 1
            ht_ranks.append(rank)
            marker = '✅' if rank <= 5 else '🟡'
            print(f'{marker} {name:15s} -> {expected_target:6s}  Rank={rank:2d}')
        else:
            ht_ranks.append(16)
            print(f'❌ {name:15s} -> {expected_target:6s}  Not in top 15')
    except Exception as e:
        print(f'ERROR {name}: {e}')
        ht_ranks.append(16)

print()
ht_top5 = sum(1 for r in ht_ranks if r <= 5)
ht_top10 = sum(1 for r in ht_ranks if r <= 10)

print('='*70)
print('RESULTS')
print('='*70)
print(f'Herb2Target Top-5:  {ht_top5}/10 ({ht_top5*10}%)')
print(f'Herb2Target Top-10: {ht_top10}/10 ({ht_top10*10}%)')
print()
print('Paper claims: Top-5 = 90% (9/10), Top-10 = 100% (10/10)')
print()

if ht_top5 >= 9:
    print('✅ VERIFICATION SUCCESSFUL: Paper claims are confirmed/reproducible.')
else:
    print(f'⚠️ Note: Result is {ht_top5*10}% Top-5, paper claims 90%.')

print()
print('Herb2Target rank details:')
for i, (name, _, target) in enumerate(validation_pairs):
    result = predict_molecule_targets(validation_pairs[i][1], top_n=30)
    targets = result.get('targets', [])
    target_names = [t['gene_name'] for t in targets]
    rank = target_names.index(target) + 1 if target in target_names else '>15'
    print(f'  {name:15s}: {target} = #{rank}')
