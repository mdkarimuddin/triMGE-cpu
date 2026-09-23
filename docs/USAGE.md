# Usage (individual CLIs)

Prefer the one-shot workflow in [EASY_PIPELINE.md](EASY_PIPELINE.md) when possible.  
This page documents each command separately.

## 1. Classify — `trimge-cpu`

```bash
trimge-cpu -i sequences.fasta -o predictions.tsv --models-dir models
```

JSON output:

```bash
trimge-cpu -i sequences.fasta -o predictions.json --format json --models-dir models
```

| Option | Default | Description |
|--------|---------|-------------|
| `-i` / `--input` | required | Input FASTA |
| `-o` / `--output` | required | Output TSV or JSON |
| `--format` | `tsv` | `tsv` or `json` |
| `--models-dir` | `./models` | Exported model directory |
| `--verbose` | off | Debug logging |

Equivalent:

```bash
python -m trimgecpu -i sequences.fasta -o predictions.tsv --models-dir models
```

### Two-stage logic

1. Compute k-mer k5 features for every sequence  
2. SVM → `phage` / `phage_plasmid` / `plasmid`  
3. If `phage_plasmid` → second SVM → one of 15 subtypes  

## 2. Combine FASTAs — `trimge-cpu-combine`

```bash
trimge-cpu-combine --input-dir /path/to/fastas -o combined.fna
```

Or:

```bash
trimge-cpu-combine -i a.fna b.fna -o combined.fna
```

| Option | Description |
|--------|-------------|
| `--recursive` | Recurse under `--input-dir` |
| `--no-prefix` | Do not prefix IDs with filename |
| `--separator` | Prefix separator (default `\|`) |
| `--min-length` | Drop short sequences (default 100) |
| `--manifest` | Manifest TSV path |

By default headers become `filename_stem|original_id` to avoid collisions.

## 3. Reports — `trimge-cpu-report`

Requires a predictions TSV from `trimge-cpu`:

```bash
trimge-cpu-report -i predictions.tsv -o report_out
```

| Option | Default | Description |
|--------|---------|-------------|
| `-i` / `--input` | required | Predictions TSV |
| `-o` / `--output-dir` | required | Report root directory |
| `--steps` | `1,2` | `1` = 3-class; `2` = P–P subtype |
| `--verbose` | off | Debug logging |

Examples:

```bash
# 3-class report only
trimge-cpu-report -i predictions.tsv -o report_out --steps 1

# P–P subtypes only
trimge-cpu-report -i predictions.tsv -o report_out --steps 2
```

See [REPORTS.md](REPORTS.md) for output files.

## 4. Full pipeline — `trimge-cpu-pipeline`

See [EASY_PIPELINE.md](EASY_PIPELINE.md).
