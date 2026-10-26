#!/usr/bin/env python3
"""Submit the isolated 10k triplet array and its dependent correction job."""
from datetime import datetime, timezone
import csv
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess
import os
import sys

ROOT = Path(__file__).resolve().parent


def sbatch(*args):
    result = subprocess.run(["sbatch", "--parsable", *map(str, args)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True, check=True)
    job_id = result.stdout.strip().split(";")[0]
    if not job_id.isdigit():
        raise ValueError(f"Unexpected sbatch output: {result.stdout}")
    return job_id


def main():
    with (ROOT / ".submission.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        release = next(p for p in ROOT.parents if (p / '.epigenetic_oinfo_root').is_file())
        os.environ['EPINFO_ROOT'] = str(release)
        # Verify the current shared implementation, not the obsolete snapshot file list.
        subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', str(release / 'tests'), '-v'], check=True)
        record_path = ROOT / "submission.json"
        record = json.loads(record_path.read_text()) if record_path.exists() else dict(
            n_samples=82, n_permutations=10000, permutation_min_order=3,
            permutation_max_order=3, cpus_per_task=32, max_parallel_samples=None,
            time_limit="24:00:00",
            output_root=str(ROOT))
        if "array_job_id" not in record:
            record["array_job_id"] = sbatch(ROOT / "submit_array.sh")
            record["array_submitted_at_utc"] = datetime.now(timezone.utc).isoformat()
            record_path.write_text(json.dumps(record, indent=2) + "\n")
            print("Submitted array:", record["array_job_id"], flush=True)
        else:
            print("Array already submitted:", record["array_job_id"], flush=True)
        with (ROOT / "input_manifest.tsv").open() as handle:
            samples = list(csv.DictReader(handle, delimiter="\t"))
        with (ROOT / "task_submissions.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "array_task_id", "task_index", "output_dir"])
            for entry in samples:
                index = entry["task_index"]
                writer.writerow([entry["sample"], record["array_job_id"] + "_" + index,
                                 index, entry["output_dir"]])
        if "correction_job_id" not in record:
            record["correction_job_id"] = sbatch(
                "--dependency=afterok:" + record["array_job_id"], ROOT / "submit_correction.sh")
            record["correction_submitted_at_utc"] = datetime.now(timezone.utc).isoformat()
            record_path.write_text(json.dumps(record, indent=2) + "\n")
            print("Submitted dependent correction:", record["correction_job_id"], flush=True)
        else:
            print("Correction already submitted:", record["correction_job_id"], flush=True)


if __name__ == "__main__":
    main()
