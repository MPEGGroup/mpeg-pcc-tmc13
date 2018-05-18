#!/bin/bash
#
# Generate a configuration tree in $PWD from YAML files in the same
# directory.

set -e

script_dir=$(dirname $0)

for f in ctc-*.yaml
do
	$script_dir/gen-cfg.pl --no-skip-sequences-without-src $f
	rm config-merged.yaml
done
