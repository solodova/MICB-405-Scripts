#!/bin/bash

echo "Running master script"

echo "Running bwa alignment"
./bwa.sh
echo "bwa alignment complete"

echo "Running sambamba tools"
./sambamba.sh
echo "sambamba tools complete"

echo "Running macs2 peaks"
./macs2_peaks.sh
echo "macs2 peaks complete"

echo "Running macs2 bdgdiff"
./macs2_diff.sh
echo "macs2 bdgdiff complete"

echo "master script complete"