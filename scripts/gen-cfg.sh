#!/bin/bash
#
# Generate a configuration tree in $PWD from YAML files in the same
# directory.

set -e
shopt -s nullglob

script_dir=$(dirname $0)
src_cfg_dir="${1}${1:+/}"
shift || true;

for f in ${src_cfg_dir}cat3-tmc13-ctc-*.yaml
do
	echo $f ...
	$script_dir/gen-cfg.pl --no-skip-sequences-without-src "$@" $f
	rm config-merged.yaml
done


# NB, cat1 configs need the sequence parameters
for f in ${src_cfg_dir}cat1-tmc13-ctc-*.yaml
do
	echo $f ...
	$script_dir/gen-cfg.pl --no-skip-sequences-without-src "$@" $f ${src_cfg_dir}sequences-cat1.yaml
	rm config-merged.yaml
done
