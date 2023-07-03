#!/bin/bash
#
# Generate a configuration tree in $PWD from YAML files in the same
# directory.

set -e
shopt -s nullglob

script_dir="$(dirname $0)"
geom="octree"
attr="predlift"
pred="intra"
src_cfg_dir=""

while (( $# )); do
	case $1 in
	--octree) geom="octree" ;;
	--predgeom) geom="predgeom" ;;
	--trisoup) geom="trisoup" ;;
	--raht) attr="raht" ;;
	--pred-lift) attr="predlift" ;;
	--intra) pred="intra" ;;
	--inter) pred="inter" ;;
	--all) all=1 ;;
	--cfgdir=*) src_cfg_dir="${1#--cfgdir=}/" ;;
	--) shift; break ;;
	--help|*)
		echo -e "usage:\n $0\n" \
			"    [[--octree|--predgeom|--trisoup] [--raht|--pred-lift] [--intra|--inter]  | --all]\n" \
			"    [--cfgdir=<dir>]"
		exit 1
	esac
	shift;
done

extra_args=("$@")

inter_sequences=(
	ford_01_q1mm
	ford_02_q1mm
	ford_03_q1mm
	qnxadas-junction-approach
	qnxadas-junction-exit
	qnxadas-motorway-join
	qnxadas-navigating-bends
	innovizQC1
	innovizQC2
	innovizQC3
)

# dirty hard mapping of motion files to sequences
declare -A inter_sequences_motion_octree=(
	[ford_01_q1mm]=../global-motion-files/globalMotion/ford_01_q1mm-global-motion-matrix-estimated.txt
	[ford_02_q1mm]=../global-motion-files/globalMotion/ford_02_q1mm-global-motion-matrix-estimated.txt
	[ford_03_q1mm]=../global-motion-files/globalMotion/ford_03_q1mm-global-motion-matrix-estimated.txt
	[qnxadas-junction-approach]=../global-motion-files/globalMotion/qnxadas-junction-approach-global-motion-matrix-estimated.txt
	[qnxadas-junction-exit]=../global-motion-files/globalMotion/qnxadas-junction-exit-global-motion-matrix-estimated.txt
	[qnxadas-motorway-join]=../global-motion-files/globalMotion/qnxadas-motorway-join-global-motion-matrix-estimated.txt
	[qnxadas-navigating-bends]=../global-motion-files/globalMotion/qnxadas-navigating-bends-global-motion-matrix-estimated.txt
	[innovizQC1]=../global-motion-files/globalMotion/gps-innovizqc1-matrix.txt
	[innovizQC2]=../global-motion-files/globalMotion/gps-innovizqc2-matrix.txt
	[innovizQC3]=../global-motion-files/globalMotion/gps-innovizqc3-matrix.txt
)
declare -A inter_sequences_motion_predgeom=(
	[ford_01_q1mm]=../global-motion-files/globalMotion/ford_01_q1mm-global-motion-matrix-estimated.txt
	[ford_02_q1mm]=../global-motion-files/globalMotion/ford_02_q1mm-global-motion-matrix-estimated.txt
	[ford_03_q1mm]=../global-motion-files/globalMotion/ford_03_q1mm-global-motion-matrix-estimated.txt
	[qnxadas-junction-approach]=../global-motion-files/globalMotion/qnxadas-junction-approach-global-motion-matrix-estimated.txt
	[qnxadas-junction-exit]=../global-motion-files/globalMotion/qnxadas-junction-exit-global-motion-matrix-estimated.txt
	[qnxadas-motorway-join]=../global-motion-files/globalMotion/qnxadas-motorway-join-global-motion-matrix-estimated.txt
	[qnxadas-navigating-bends]=../global-motion-files/globalMotion/qnxadas-navigating-bends-global-motion-matrix-estimated.txt
)

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
	octree-raht-ctc-lossless-geom-lossless-attrs.yaml
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
	local pred=$3

	old_src_cfg_dir=${src_cfg_dir}
	seq_src_cfg_dir=${src_cfg_dir}
	if [[ $pred != "inter" ]]
	then
		outdir="$geom-$attr/"
	else
		outdir="$geom-$attr-inter/"
		src_cfg_dir=${src_cfg_dir}inter/
	fi
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
		"${seq_src_cfg_dir}sequences-cat1.yaml" \
		"${seq_src_cfg_dir}sequences-cat3.yaml" \
		"${extra_args[@]}"

	rm -f "$outdir/config-merged.yaml"
	src_cfg_dir=${old_src_cfg_dir}

	if [[ "$pred" == "inter" ]]
	then
		for seq in ${inter_sequences[@]}
		do
			for d in $(find $outdir -name $seq -type d)
			do
				local vars=$(find $d/* -type d)
				local item="inter_sequences_motion_${geom}[$seq]"
				local src_motion="${!item}"
				if [[ "$vars" == "" ]]
				then
					[ -z "$src_motion" ] || cp ${src_cfg_dir}$src_motion $d
				else
					for var in $vars
					do
						[ -z "$src_motion" ] || cp ${src_cfg_dir}$src_motion $var
					done
				fi
			done
		done
	fi
}

if [[ "$all" != "1" ]]
then
	do_one_cfgset "$geom" "$attr" "$pred"
else
	do_one_cfgset "octree" "predlift" "intra"
	do_one_cfgset "octree" "raht" "intra"
	do_one_cfgset "predgeom" "predlift" "intra"
	do_one_cfgset "predgeom" "raht" "intra"
	do_one_cfgset "trisoup" "predlift" "intra"
	do_one_cfgset "trisoup" "raht" "intra"
	do_one_cfgset "octree" "predlift" "inter"
	do_one_cfgset "octree" "raht" "inter"
	do_one_cfgset "predgeom" "predlift" "inter"
	do_one_cfgset "predgeom" "raht" "inter"
fi
