import pandas as pd
from pathlib import Path
import shutil
import warnings
warnings.filterwarnings("ignore")

# Import your prediction function from the same folder
from .predict_psr import predict_psr_on_clone

def run_ml_prediction(folder: Path):
    folder = Path(folder)

    # Find only uncompressed clone files
    clone_files = list(folder.glob("*clones.csv")) + list(folder.glob("*clones.xlsx"))

    if not clone_files:
        print("No uncompressed clone files found — Step 6 skipped.")
        return

    print(f"Found {len(clone_files)} clone files for PSR prediction")

    for f in clone_files:
        target = f.stem.replace("_clone", "")
        print(f"\nProcessing target: {target} ({f.suffix})")
       
        # Make backup
        backup_file = folder / f"{target}_clone_backup{f.suffix}"
        if not backup_file.exists():
            shutil.copy2(f, backup_file)
            print(f"  Backup created: {backup_file.name}")

        # Run prediction (modifies file in place)
        try:
        
            predict_psr_on_clone(str(f))
            print(f"  PSR prediction completed on {f.name}")
        except Exception as e:
            print(f"  Prediction failed for {f.name}: {e}")
            continue

        print(f"  Updated original file with PSR predictions: {f.name}")

    print("\nStep 6 complete! PSR predictions added directly to clone files")
    print("Backups created as *_clone_backup.csv / *_clone_backup.xlsx")
