ply-merge: A tool to merge/split ply frames
===========================================

The ply-merge tool combines point clouds from multiple ply files into a
single output with an extra per-attribute frameindex property that
identifies which input frame each point belongs to.  The tool is also
able to reverse the process and split a merged point cloud into
individual frames.

Merge operation
---------------
From a sequence of input ply files, the merge mode repeatedly takes
`groupSize` frames, merges their contents and produces a single output
file with an additional frameindex property.

Split operation
---------------
From a sequence of input ply files, and for each value of the frameindex
property, the split mode extracts all points with the same frameindex.


Options
-------

### `--help`
Print a list of available command line (and configuration file) options
along with their default values and exit.


### `--config=FILE`, `-c`
Specifies a configuration file to be immediately loaded.


### `--mode=MODE`
Selects the mode of operation.

  | Value | Description                                |
  |:-----:| ------------------------------------------ |
  | merge | Combines multiple input files into outputs |
  | split | Splits input files into multiple outputs   |


### `--outputBinaryPly=0|1`
Sets the output format of PLY files (Binary=1, ASCII=0).  Reading and
writing binary PLY files is more efficient than the ASCII variant,
but are less suited to simple scripts and direct human inspection.

If outputting non-integer point co-ordinates (eg, due to the output
geometry scaling), the precision of the binary and ASCII versions are
not identical.


### `--srcPath=FILESPEC`, `--outPath=FILESPEC`
Specify the input and output ply file names.  The tool respectively
replaces any instance of a '%d' printf format directive with the
current input/output frame numbers.


### `--firstFrameNum=INT-VALUE`, `--firstOutputFrameNum=INT-VALUE`
The initial frame number of the respective input/output sequences.


### `--frameCount=INT-VALUE`
The number of input frames to process


### `--groupSize=INT-VALUE`
(Merge mode only)
The number of input frames to merge into each output frame.


Example usage
-------------

### Merging

```console
$ build/tmc3/ply-merge --mode=merge \
    --srcPath=path/to/Ford_01_q_1mm/Ford_01_vox1mm-%.04d.ply \
    --outPath=ford-01-vox1mm-merge8f-%.04d.ply \
    --firstFrameNum=100 \
    --frameCount=32
MPEG PCC ply merge/split tool from Test Model C13
    help                : 0
    config              : ...
    mode                : merge
    srcPath             : "path/to/Ford_01_q_1mm/Ford_01_vox1mm-%.04d.ply"
    outPath             : "ford-01-vox1mm-merge8f-%.04d.ply"
    outputBinaryPly     : 0
    firstFrameNum       : 100
    firstOutputFrameNum : 0
    frameCount          : 32
    groupSize           : 8
ford-01-vox1mm-merge8f-0000.ply
ford-01-vox1mm-merge8f-0001.ply
ford-01-vox1mm-merge8f-0002.ply
ford-01-vox1mm-merge8f-0003.ply
```

### Splitting

```console
$ build/tmc3/ply-merge --mode=split \
    --srcPath=ford-01-vox1mm-merge8f-%.04d.ply \
    --outPath=split-Ford_01_vox1mm-%.04d.ply \
    --firstOutputFrameNum=100 \
    --frameCount=4
```
    help                : 0
    config              : ...
    mode                : split
    srcPath             : "ford-01-vox1mm-merge8f-%.04d.ply"
    outPath             : "split-Ford_01_vox1mm-%.04d.ply"
    outputBinaryPly     : 0
    firstFrameNum       : 0
    firstOutputFrameNum : 100
    frameCount          : 32
    groupSize           : 8
split-Ford_01_vox1mm-0100.ply
split-Ford_01_vox1mm-0101.ply
split-Ford_01_vox1mm-0102.ply
...
split-Ford_01_vox1mm-0131.ply
```
