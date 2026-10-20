from __future__ import annotations

import re
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .common import chrom_key, open_text, read_json, read_list, safe_name, write_json


CONTIG_RE = re.compile(r"##contig=<ID=([^,>]+).*?(?:length|Length)=(\d+)")


def _parse_gt(sample: str, gt_index: int) -> int:
    if gt_index == 0:
        gt = sample.split(":", 1)[0]
    else:
        parts = sample.split(":")
        if gt_index >= len(parts):
            return -1
        gt = parts[gt_index]
    if len(gt) >= 3 and gt[1] in "/|":
        a, b = gt[0], gt[2]
        if a in "01" and b in "01":
            return (a == "1") + (b == "1")
    if gt == "0":
        return 0
    if gt == "1":
        return 2
    return -1


def _source_fingerprint(path):
    source = Path(path).resolve()
    stat = source.stat()
    return {"path": str(source), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def build_cache(vcf, cultivar_list, wild_list, outdir, resume=False, force=False, chunk_sites=10000):
    out = Path(outdir).resolve()
    manifest_path = out / "manifest.json"
    if resume and manifest_path.exists():
        manifest = read_json(manifest_path)
        cultivars = read_list(cultivar_list)
        wild = read_list(wild_list) if wild_list else []
        if manifest.get("source_fingerprint") != _source_fingerprint(vcf):
            raise SystemExit("Cached VCF fingerprint differs from --vcf; use a new --cache-dir or --force")
        if manifest.get("cultivar_samples") != cultivars or manifest.get("wild_samples") != wild:
            raise SystemExit("Cached sample lists differ from current lists; use a new --cache-dir or --force")
        return manifest
    if out.exists():
        if not force:
            raise SystemExit(f"Cache directory already exists: {out}. Use --resume to validate/reuse it or --force to rebuild it")
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)

    cultivars = read_list(cultivar_list)
    wild = read_list(wild_list) if wild_list else []
    if not cultivars:
        raise SystemExit("cultivar list is empty")
    overlap = sorted(set(cultivars) & set(wild))
    if overlap:
        raise SystemExit("Samples occur in both cultivar and wild lists: " + ",".join(overlap[:20]))
    wanted = cultivars + wild
    called_counts = np.zeros(len(wanted), dtype=np.int64)
    missing_counts = np.zeros(len(wanted), dtype=np.int64)
    heterozygous_counts = np.zeros(len(wanted), dtype=np.int64)
    contig_lengths = {}
    chroms = {}
    selected = None
    current_chrom = None
    current_meta = None
    handles = {}
    positions, genotypes, mafs, acs, ans, refalt = [], [], [], [], [], bytearray()
    completed_chroms = set()
    last_pos = -1
    n_input = n_kept = 0
    started = time.time()

    def close_handles():
        nonlocal handles
        for handle in handles.values():
            handle.close()
        handles = {}

    def flush():
        nonlocal positions, genotypes, mafs, acs, ans, refalt
        if not positions:
            return
        np.asarray(positions, dtype=np.int64).tofile(handles["positions"])
        np.asarray(genotypes, dtype=np.int8).tofile(handles["genotypes"])
        np.asarray(mafs, dtype=np.float32).tofile(handles["maf"])
        np.asarray(acs, dtype=np.int32).tofile(handles["ac"])
        np.asarray(ans, dtype=np.int32).tofile(handles["an"])
        handles["refalt"].write(refalt)
        positions, genotypes, mafs, acs, ans, refalt = [], [], [], [], [], bytearray()

    def start_chrom(chrom):
        nonlocal handles, current_meta
        safe = safe_name(chrom)
        cdir = out / safe
        cdir.mkdir(parents=True, exist_ok=True)
        handles = {
            "positions": (cdir / "positions.i64").open("ab"),
            "genotypes": (cdir / "genotypes.i8").open("ab"),
            "maf": (cdir / "cultivar_maf.f32").open("ab"),
            "ac": (cdir / "cultivar_ac.i32").open("ab"),
            "an": (cdir / "cultivar_an.i32").open("ab"),
            "refalt": (cdir / "refalt.u8").open("ab"),
        }
        current_meta = {"safe": safe, "n_sites": 0, "max_pos": 0, "length": contig_lengths.get(chrom, 0)}
        chroms[chrom] = current_meta

    with open_text(vcf) as handle:
        for line in handle:
            if line.startswith("##contig=<"):
                match = CONTIG_RE.search(line)
                if match:
                    contig_lengths[match.group(1)] = int(match.group(2))
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM\t"):
                header_samples = line.rstrip("\r\n").split("\t")[9:]
                lookup = {sample: i for i, sample in enumerate(header_samples)}
                missing = [sample for sample in wanted if sample not in lookup]
                if missing:
                    raise SystemExit("Requested samples absent from VCF: " + ",".join(missing[:30]))
                selected = [lookup[sample] + 9 for sample in wanted]
                continue
            if line.startswith("#"):
                continue
            if selected is None:
                raise SystemExit("VCF data encountered before #CHROM header")
            n_input += 1
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 10 or len(fields[3]) != 1 or len(fields[4]) != 1 or fields[4] == "." or "," in fields[4]:
                continue
            chrom = fields[0]
            pos = int(fields[1])
            if chrom != current_chrom:
                if current_chrom is not None:
                    flush(); close_handles(); completed_chroms.add(current_chrom)
                if chrom in completed_chroms:
                    raise SystemExit(f"VCF is not chromosome-sorted: {chrom} occurs in multiple runs")
                current_chrom = chrom
                last_pos = -1
                start_chrom(chrom)
            if pos < last_pos:
                raise SystemExit(f"VCF is not position-sorted within {chrom}: {pos} after {last_pos}")
            last_pos = pos
            fmt = fields[8].split(":")
            if "GT" not in fmt:
                continue
            gt_index = fmt.index("GT")
            row = np.fromiter((_parse_gt(fields[column], gt_index) for column in selected), dtype=np.int8, count=len(selected))
            called_counts += row >= 0
            missing_counts += row < 0
            heterozygous_counts += row == 1
            cultivar_gt = row[:len(cultivars)]
            called = cultivar_gt >= 0
            an = int(called.sum()) * 2
            ac = int(cultivar_gt[called].sum()) if an else 0
            af = ac / an if an else np.nan
            maf = min(af, 1.0 - af) if an else np.nan
            positions.append(pos - 1)
            genotypes.append(row)
            mafs.append(maf)
            acs.append(ac)
            ans.append(an)
            refalt.extend(fields[3].encode("ascii", "replace")[:1] + fields[4].encode("ascii", "replace")[:1])
            current_meta["n_sites"] += 1
            current_meta["max_pos"] = pos
            n_kept += 1
            if len(positions) >= chunk_sites:
                flush()
            if n_input % 500000 == 0:
                print(f"[cache] input={n_input:,} kept={n_kept:,} elapsed={time.time()-started:.1f}s", file=sys.stderr, flush=True)
    if current_chrom is not None:
        flush(); close_handles()
    if not chroms:
        raise SystemExit("No biallelic SNP records were retained")
    for chrom, meta in chroms.items():
        meta["length"] = max(meta["length"], meta["max_pos"])
    manifest = {
        "version": 1,
        "source_vcf": str(Path(vcf).resolve()),
        "source_fingerprint": _source_fingerprint(vcf),
        "cultivar_samples": cultivars,
        "wild_samples": wild,
        "samples": wanted,
        "n_cultivar": len(cultivars),
        "n_wild": len(wild),
        "n_input_records": n_input,
        "n_retained_biallelic_snps": n_kept,
        "sample_qc": [
            {
                "sample": sample,
                "group": "CULTIVAR" if index < len(cultivars) else "WILD",
                "called_genotypes": int(called_counts[index]),
                "missing_genotypes": int(missing_counts[index]),
                "call_rate": float(called_counts[index] / max(1, n_kept)),
                "heterozygous_genotypes": int(heterozygous_counts[index]),
                "heterozygosity_among_called": float(heterozygous_counts[index] / max(1, called_counts[index])),
            }
            for index, sample in enumerate(wanted)
        ],
        "chromosomes": {chrom: chroms[chrom] for chrom in sorted(chroms, key=chrom_key)},
    }
    write_json(manifest_path, manifest)
    return manifest


@dataclass
class ChromCache:
    root: Path
    chrom: str
    meta: dict
    n_samples: int
    n_cultivar: int

    def __post_init__(self):
        directory = self.root / self.meta["safe"]
        n = int(self.meta["n_sites"])
        self.positions = np.memmap(directory / "positions.i64", dtype=np.int64, mode="r", shape=(n,))
        self.genotypes = np.memmap(directory / "genotypes.i8", dtype=np.int8, mode="r", shape=(n, self.n_samples))
        self.maf = np.memmap(directory / "cultivar_maf.f32", dtype=np.float32, mode="r", shape=(n,))
        self.ac = np.memmap(directory / "cultivar_ac.i32", dtype=np.int32, mode="r", shape=(n,))
        self.an = np.memmap(directory / "cultivar_an.i32", dtype=np.int32, mode="r", shape=(n,))
        self.refalt = np.memmap(directory / "refalt.u8", dtype=np.uint8, mode="r", shape=(n, 2))

    def bounds(self, start0: int, end0: int):
        return int(np.searchsorted(self.positions, start0, side="left")), int(np.searchsorted(self.positions, end0, side="left"))


def open_cache(cache_dir):
    root = Path(cache_dir).resolve()
    manifest = read_json(root / "manifest.json")
    caches = {
        chrom: ChromCache(root, chrom, meta, len(manifest["samples"]), manifest["n_cultivar"])
        for chrom, meta in manifest["chromosomes"].items()
    }
    return manifest, caches
