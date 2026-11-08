#!/usr/bin/env python3
"""Run a one-design RFantibody installation smoke test in Linux/WSL.

From the repository root: python3 src/generate/pilot_rfantibody.py
Uses chain T (lysozyme) from the committed 1MEL test fixture and the upstream
NbBCII10 framework. This is an execution check, not a binding-quality benchmark.
"""
import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time


ROOT = Path(__file__).resolve().parents[2]


def inspect_pdb(path):
    """Check finite coordinates, H/T residues, CDR labels, and optional scores."""
    residues = {}
    labels = set()
    scores = {}
    for line in path.read_text().splitlines():
        if line.startswith("ATOM"):
            xyz = [float(line[start:start + 8]) for start in (30, 38, 46)]
            if not all(math.isfinite(x) for x in xyz):
                raise ValueError(f"Nonfinite coordinates in {path}")
            if line[12:16].strip() == "CA":
                key = (line[21], line[22:27])
                if key in residues:
                    raise ValueError(f"Duplicate CA residue in {path}: {key}")
                residues[key] = line[17:20]
        elif line.startswith("REMARK PDBinfo-LABEL:"):
            labels.update(line.split()[3:])
        elif line.startswith("SCORE "):
            key, value = line[6:].split(":", 1)
            scores[key] = float(value)
    chains = {chain: sum(k[0] == chain for k in residues)
              for chain in {k[0] for k in residues}}
    if not chains.get("H") or not chains.get("T"):
        raise ValueError(f"Missing H or T residues in {path}: {chains}")
    if not {"H1", "H2", "H3"}.issubset(labels):
        raise ValueError(f"Missing CDR labels in {path}: {labels}")
    for key, value in scores.items():
        if not math.isfinite(value):
            if (not chains.get("L") and math.isnan(value) and key in
                    {f"framework_aligned_L{i}_rmsd" for i in (1, 2, 3)}):
                scores[key] = None  # RF2 computes RMSD on empty light-chain masks.
            else:
                raise ValueError(f"Nonfinite score {key} in {path}")
    return {"file": str(path), "residues": chains, "scores": scores}


def validate_outputs(out, summary):
    for key, pattern in (("backbones", "*.pdb"), ("sequences", "*.pdb"),
                         ("predictions", "*best.pdb")):
        files = sorted((out / key).glob(pattern))
        if len(files) != 1:
            raise ValueError(f"Expected one {key} structure, got {len(files)}")
        summary[key] = [inspect_pdb(p) for p in files]
    scores = summary["predictions"][0]["scores"]
    for name in ("interaction_pae", "pred_lddt"):
        if name not in scores or scores[name] is None:
            raise ValueError(f"Missing RF2 confidence metric: {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rfantibody-root", type=Path, default=ROOT.parent / "RFantibody")
    parser.add_argument("--output-dir", type=Path,
                        help="New or empty directory; default runs/rfantibody_pilot_<UTC timestamp>")
    parser.add_argument("--validate-only", action="store_true",
                        help="Recheck an existing --output-dir without running models")
    args = parser.parse_args()
    if args.validate_only:
        if args.output_dir is None:
            parser.error("--validate-only requires --output-dir")
        out = args.output_dir.resolve()
        summary = json.loads((out / "summary.json").read_text())
        summary.setdefault("validation_history", []).append(
            {"status": summary["status"], "error": summary.pop("error", None)})
        try:
            stages = {s["name"]: s for s in summary["stages"]}
            for name in ("gpu_check", "rfdiffusion", "proteinmpnn", "rf2"):
                if stages.get(name, {}).get("returncode") != 0:
                    raise ValueError(f"Stage {name} did not complete successfully")
            validate_outputs(out, summary)
            summary["status"] = "passed"
        except Exception as exc:
            summary["status"] = "failed"
            summary["error"] = f"{type(exc).__name__}: {exc}"
        summary["validated_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
        result = f"Revalidation: {summary['status']}. {summary.get('error', '')}"
        with (out / "pilot.log").open("a") as log:
            log.write(result + "\n")
        print(result)
        return 0 if summary["status"] == "passed" else 1
    upstream = args.rfantibody_root.resolve()
    out = (args.output_dir or ROOT / "runs" / (
        "rfantibody_pilot_" + dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))).resolve()
    if out.exists() and any(out.iterdir()):
        parser.error(f"Output directory must be empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    summary = {"status": "running", "purpose": "installation smoke test, not binding validation",
               "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
               "rfantibody_root": str(upstream), "output_dir": str(out), "stages": [],
               "parameters": {"backbones": 1, "sequences": 1, "diffusion_steps": 50,
                              "loops": "H1:7,H2:6,H3:8", "rf2_recycles": 1,
                              "temperature": 0.1, "rf2_seed": 37}}
    started = time.monotonic()
    env = os.environ.copy()
    env["PATH"] = str(upstream / ".venv/bin") + os.pathsep + env.get("PATH", "")
    env.update(RFANTIBODY_ROOT=str(upstream), PYTHONUNBUFFERED="1", HYDRA_FULL_ERROR="1")
    log_path = out / "pilot.log"

    def save():
        (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    with log_path.open("w", buffering=1) as log:
        def message(text):
            print(text, flush=True)
            log.write(text + "\n")

        def run(name, command):
            stage = {"name": name, "command": command, "status": "running"}
            summary["stages"].append(stage)
            save()
            message(f"\n[{name}] {subprocess.list2cmdline(command)}")
            tick = time.monotonic()
            try:
                # Stream output to disk: large model logs do not accumulate in RAM.
                with subprocess.Popen(command, cwd=out, env=env, stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT, text=True, errors="replace") as proc:
                    for line in proc.stdout:
                        log.write(line)
                    code = proc.wait()
                stage["returncode"] = code
                if code:
                    raise RuntimeError(f"{name} exited with code {code}; see {log_path}")
                stage["status"] = "passed"
            except BaseException:
                stage["status"] = "failed"
                raise
            finally:
                stage["elapsed_seconds"] = round(time.monotonic() - tick, 2)
                save()
            message(f"[{name}] completed in {stage['elapsed_seconds']} seconds")

        try:
            message(f"Pilot output: {out}")
            for name in ("python", "rfdiffusion", "proteinmpnn", "rf2"):
                if not (upstream / ".venv/bin" / name).is_file():
                    raise FileNotFoundError(f"Missing installed command: {name}")
            for name in ("RFdiffusion_Ab.pt", "ProteinMPNN_v48_noise_0.2.pt", "RF2_ab.pt"):
                path = upstream / "weights" / name
                if not path.is_file() or path.stat().st_size == 0:
                    raise FileNotFoundError(f"Missing weights: {path}")
            summary["upstream_commit"] = subprocess.check_output(
                ["git", "-C", str(upstream), "rev-parse", "HEAD"], text=True).strip()
            run("gpu_check", ["python", "-c", "import torch, dgl; "
                "assert torch.cuda.is_available(), 'CUDA unavailable'; "
                "print(torch.__version__, dgl.__version__, torch.cuda.get_device_name(0)); "
                "print(torch.ones(1, device='cuda').sum().item())"])
            inputs = out / "inputs"
            inputs.mkdir()
            fixture = ROOT / "src/evaluate/testdata/hlt/1MEL.pdb"
            target = inputs / "lysozyme_T.pdb"
            atoms = [line for line in fixture.read_text().splitlines()
                     if line.startswith("ATOM") and line[21] == "T"]
            if not atoms:
                raise ValueError("No lysozyme atoms in fixture")
            target.write_text("\n".join(atoms) + "\nTER\nEND\n")
            framework = inputs / "h-NbBCII10.pdb"
            shutil.copyfile(upstream / "scripts/examples/example_inputs/h-NbBCII10.pdb", framework)
            summary["inputs"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in (target, framework)}
            backbones, sequences, predictions = [out / name for name in
                                                ("backbones", "sequences", "predictions")]
            for folder in (backbones, sequences, predictions):
                folder.mkdir()
            run("rfdiffusion", ["rfdiffusion", "-t", str(target), "-f", str(framework),
                "-o", str(backbones / "pilot"), "-n", "1", "-l", "H1:7,H2:6,H3:8",
                "--diffuser-t", "50", "--deterministic", "--no-trajectory"])
            files = sorted(backbones.glob("*.pdb"))
            if len(files) != 1:
                raise ValueError(f"Expected one backbone, got {len(files)}")
            summary["backbones"] = [inspect_pdb(p) for p in files]
            run("proteinmpnn", ["proteinmpnn", "-i", str(backbones), "-o", str(sequences),
                "-n", "1", "-l", "H1,H2,H3", "-t", "0.1", "--deterministic", "--debug"])
            files = sorted(sequences.glob("*.pdb"))
            if len(files) != 1:
                raise ValueError(f"Expected one sequence structure, got {len(files)}")
            summary["sequences"] = [inspect_pdb(p) for p in files]
            run("rf2", ["rf2", "-i", str(sequences), "-o", str(predictions),
                "-r", "1", "-s", "37", "--hotspot-show-prop", "0.0"])
            validate_outputs(out, summary)
            summary["status"] = "passed"
            message("PASS: all three stages produced structures with H/T chains, CDR labels, "
                    "finite coordinates, and RF2 confidence metrics.")
        except (Exception, KeyboardInterrupt) as exc:
            summary["status"] = "failed"
            summary["error"] = f"{type(exc).__name__}: {exc}"
            message("FAIL: " + summary["error"])
        finally:
            summary["elapsed_seconds"] = round(time.monotonic() - started, 2)
            save()
            message(f"Summary: {out / 'summary.json'}")
    return 0 if summary["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
