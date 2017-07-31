#!/usr/bin/env python3
import yaml
import glob
from os import path
import sys

if len(sys.argv) != 2:
    print("USAGE: aggregate_qc REPORT_DIR")
    exit(1)

repdir = sys.argv[1]
print("name", "nread", "adapt_trimmed", "merged", "qc_trimmed", "qc_dropped", sep='\t')
for sample in glob.glob(path.join(repdir, "*.yml")):
    sn = path.splitext(path.basename(sample))[0]
    with open(sample) as fh:
        report = yaml.load(fh)
    nread = report[2]['AdaptorTrimPE']['output']['num_reads']
    trimmed = report[2]['AdaptorTrimPE']['output']['num_trimmed']
    merged = report[2]['AdaptorTrimPE']['output']['num_merged']
    qctrimmed = report[3]['WindowedQualTrim']['output']['num_trimmed']
    qcfilt = report[3]['WindowedQualTrim']['output']['num_dropped']
    print(sn, nread, trimmed, merged, qctrimmed, qcfilt, sep='\t')
print("Note: Numbers here are reads, not read-pairs", file=sys.stderr)
