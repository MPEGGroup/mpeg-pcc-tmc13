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
	--predgeom) geom="predgeom" ;;
	--trisoup) geom="trisoup" ;;
	--raht) attr="raht" ;;
	--pred-lift) attr="predlift" ;;
	--all) all=1 ;;
	--cfgdir=*) src_cfg_dir="${1#--cfgdir=}/" ;;
	--) shift; break ;;
	--help|*)
		echo -e "usage:\n $0\n" \
			"    [[--octree|--predgeom|--trisoup] [--raht|--pred-lift] | --all]\n" \
			"    [--cfgdir=<dir>]"
		exit 1
	esac
	shift;
done

extra_args=("$@")

##
# NB: it is important that the configs in each config set are
# capable of being merged together by gen-cfg.pl.  Ie, no two
# configs may have different definitions of one category.
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

cfg_predgeom_predlift=(
	"${cfg_octree_predlift[@]}"
	cfg-predgeom.yaml
)

cfg_predgeom_raht=(
	"${cfg_octree_raht[@]}"
	cfg-predgeom.yaml
)

cfg_trisoup_predlift=(
	trisoup-liftt-ctc-lossy-geom-lossy-attrs.yaml
)

cfg_trisoup_raht=(
	trisoup-raht-ctc-lossy-geom-lossy-attrs.yaml
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
	done

	# NB: specifying extra_args at the end does not affect option
	# processing since gen-cfg.pl is flexible in argument positions
	$script_dir/gen-cfg.pl \
		--prefix="$outdir" --no-skip-sequences-without-src \
		"${!cfgset/#/${src_cfg_dir}}" \
		"${src_cfg_dir}sequences-cat1.yaml" \
		"${src_cfg_dir}sequences-cat3.yaml" \
		"${extra_args[@]}"

	rm -f "$outdir/config-merged.yaml"
}

if [[ "$all" != "1" ]]
then
	do_one_cfgset "$geom" "$attr"
else
	do_one_cfgset "octree" "predlift"
	do_one_cfgset "octree" "raht"
	do_one_cfgset "predgeom" "predlift"
	do_one_cfgset "predgeom" "raht"
	do_one_cfgset "trisoup" "predlift"
	do_one_cfgset "trisoup" "raht"
fi
