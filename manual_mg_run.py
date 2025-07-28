#!/usr/bin/env python3
"""
Manual MadGraph workflow - replicating what MadMiner does internally
This script shows the step-by-step process to run MadGraph manually
"""

print("""
=================================================================
MANUAL MADGRAPH WORKFLOW - Step by Step
=================================================================

This replicates what MadMiner does internally but manually.

STEP 1: Generate the Process
----------------------------
""")

mg_commands = [
    "import model SMEFTatNLO-LO",
    "generate p p > H H QED=2 QCD=2 NP=2 [QCD]", 
    "output signal_sm_manual",
    "quit"
]

print("MadGraph commands to run:")
for cmd in mg_commands:
    print(f"MG5_aMC> {cmd}")

print("""

STEP 2: Copy Cards to Generated Process
---------------------------------------
After process generation, copy your cards:
""")

copy_commands = [
    "cp cards/run_cards/run_card_signal_14TeV.dat signal_sm_manual/Cards/run_card.dat",
    "cp cards/madspin_card.dat signal_sm_manual/Cards/",
    "cp cards/pythia8_card.dat signal_sm_manual/Cards/", 
    "cp cards/restrict_default.dat signal_sm_manual/Cards/param_card.dat"
]

for cmd in copy_commands:
    print(f"$ {cmd}")

print("""

STEP 3: Launch Event Generation
-------------------------------
cd signal_sm_manual
./bin/generate_events auto

STEP 4: Interactive Options
---------------------------
When prompted, select:
- 0 or 'auto' to proceed with current settings
- This will use your NLO + PYTHIA8 + MadSpin setup

STEP 5: Monitor Output
---------------------
Watch for:
✓ Compilation successful
✓ NLO matrix element generation  
✓ PYTHIA8 shower
✓ MadSpin decays (h > b b~, h > a a)
✓ Final LHE/HEP files generated

=================================================================
""") 