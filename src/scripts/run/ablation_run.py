import os
import subprocess
import time
import re

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
MAKEFILE_PATH = os.path.join(REPO_ROOT, "Makefile")
LOG_PATH = os.path.join(REPO_ROOT, "results", "experiment_result.log")
RUN_ALL_PS1 = os.path.join(REPO_ROOT, "scripts", "run", "run_all.ps1")

combinations = [
    ("Base (Phase 5 Default)", ""),
    ("Phase 4 (AVX2 Dense)", "-DWJ_USE_AVX2_DENSE=1"),
    ("Phase 3 (Postings)", "-DENABLE_PHASE3_POSTINGS=1"),
    ("Phase 3 + Phase 4 (Postings + AVX2)", "-DENABLE_PHASE3_POSTINGS=1 -DWJ_USE_AVX2_DENSE=1"),
]

original_makefile = ""
with open(MAKEFILE_PATH, "r", encoding="utf-8") as f:
    original_makefile = f.read()

results = []

print("Starting Ablation Study (4 Threads)...")

for name, flags in combinations:
    print(f"\n--- Testing Combination: {name} ---")
    
    # Inject EXTRA_CXXFLAGS into Makefile
    new_makefile = original_makefile.replace("EXTRA_CXXFLAGS ?=", f"EXTRA_CXXFLAGS ?= {flags}")
    with open(MAKEFILE_PATH, "w", encoding="utf-8") as f:
        f.write(new_makefile)
    
    # Clear old log to avoid parsing previous runs
    if os.path.exists(LOG_PATH):
        os.remove(LOG_PATH)
    
    # Run test
    cmd = ["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", RUN_ALL_PS1, "-Runs", "1", "-Threads", "4"]
    print("Running tests...")
    subprocess.run(cmd, cwd=REPO_ROOT)
    
    # Parse log for test2.fasta average time
    with open(LOG_PATH, "r", encoding="utf-16") as f:
        log_content = f.read()
        
    match = re.search(r"test2\.fasta Avg Time Across Thresholds:\s*([\d\.]+)\s*seconds", log_content)
    if match:
        avg_time = float(match.group(1))
        results.append((name, avg_time))
        print(f"Result for {name}: {avg_time} s")
    else:
        print(f"Failed to parse result for {name}")
        results.append((name, None))

# Restore original Makefile
with open(MAKEFILE_PATH, "w", encoding="utf-8") as f:
    f.write(original_makefile)

print("\n\n=========================================")
print("Ablation Study Results (test2.fasta, 4 Threads)")
print("=========================================")
for name, avg_time in results:
    if avg_time is not None:
        print(f"{name:40}: {avg_time:.3f} s")
    else:
        print(f"{name:40}: FAILED")
