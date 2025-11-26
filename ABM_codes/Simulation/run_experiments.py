import os
import subprocess
import shutil

experiments = [
    {"name": "sm_only", "config": "EXPERIMENT_TYPE = 'sm_only'"},
    {"name": "mixed", "config": "EXPERIMENT_TYPE = 'mixed'"},
    {"name": "separated", "config": "EXPERIMENT_TYPE = 'separated'"},
]

for exp in experiments:
    print(f"\n=== Running experiment: {exp['name']} ===")
    
    # Modify variablevals.py
    with open("variablevals.py", "r") as f:
        lines = f.readlines()
    
    # Update experiment type
    for i, line in enumerate(lines):
        if "EXPERIMENT_TYPE" in line:
            lines[i] = f"{exp['config']}\n"
    
    with open("variablevals.py", "w") as f:
        f.writelines(lines)
    
    # Run simulation
    subprocess.run(["python", "oral_biofilm3d.py"])
    
    # Save results
    results_dir = f"results_{exp['name']}"
    if os.path.exists(results_dir):
        shutil.rmtree(results_dir)
    shutil.copytree("output", results_dir)
    
print("\nAll experiments completed!")