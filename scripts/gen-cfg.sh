#!/bin/bash
#
# Generate a configuration tree in $PWD from YAML files in the same
# directory.

set -e
shopt -s nullglob

script_dir="$(dirname $0)"
geom="octree"
attr="predlift"
src_cfg_dir=""

while (( $# )); do
	case $1 in
	--octree) geom="octree" ;;
	--trisoup) geom="trisoup" ;;
	--raht) attr="raht" ;;
	--pred-lift) attr="predlift" ;;
	--all) all=1 ;;
	--cfgdir=*) src_cfg_dir="${1#--cfgdir=}/" ;;
	--) break ;;
	--help|*)
		echo -e "usage:\n $0\n" \
			"    [[--octree|--trisoup] [--raht|--pred-lift] | --all]\n" \
			"    [--cfgdir=<dir>]"
		exit 1
	esac
	shift;
done

extra_args=("$@")

cfg_octree_predlift=(
	octree-liftt-ctc-lossless-geom-lossy-attrs.yaml
	octree-liftt-ctc-lossy-geom-lossy-attrs.yaml
	octree-predt-ctc-lossless-geom-lossless-attrs.yaml
	octree-predt-ctc-lossless-geom-nearlossless-attrs.yaml
)

cfg_octree_raht=(
	octree-raht-ctc-lossless-geom-lossy-attrs.yaml
	octree-raht-ctc-lossy-geom-lossy-attrs.yaml
)

cfg_trisoup_predlift=(
	trisoup-liftt-ctc-lossy-geom-lossy-attrs.yaml
)

cfg_trisoup_raht=(
	trisoup-liftt-ctc-lossy-geom-lossy-attrs.yaml
)

do_one_cfgset() {
	local geom=$1
	local attr=$2

	outdir="$geom-$attr/"
	mkdir -p "$outdir"

	cfgset="cfg_${geom}_${attr}[@]"

	for f in ${!cfgset}
	do
		echo "${src_cfg_dir}$f -> $outdir" ...

		# NB: specifying extra_args at the end does not affect option
		# processing since gen-cfg.pl is flexible in argument positions
		$script_dir/gen-cfg.pl \
			--prefix="$outdir" --no-skip-sequences-without-src \
			${src_cfg_dir}$f \
			${src_cfg_dir}sequences-cat1.yaml \
			${src_cfg_dir}sequences-cat3.yaml \
			"${extra_args[@]}"

		rm -f "$outdir/config-merged.yaml"
	done
}

if [[ "$all" != "1" ]]
then
	do_one_cfgset "$geom" "$attr"
else
	do_one_cfgset "octree" "predlift"
	do_one_cfgset "octree" "raht"
	do_one_cfgset "trisoup" "predlift"
	do_one_cfgset "trisoup" "raht"
fi
