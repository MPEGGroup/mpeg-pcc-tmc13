Using the codec
===============

```
./tmc3 [--help] [-c config.cfg] [--parameter=value]
```

The encoder takes as input one or more PLY files describing a point
cloud sequence with integer positions and, optionally, per-point integer
colour and reflectance attributes.

The output of the encoder is a binary bitstream encapsulated using the
G-PCC annex-B format.

Conversely, the decoder takes as input a compressed bitstream file in
G-PCC annex-B format and produces one or more reconstructed PLY file
with position and any present attribute values.

The software may be configured using either command line arguments or from
a configuration file specified using the `-c|--config=` option.

Sample configuration files are provided in the cfg/ directory.  The
utility <scripts/gen-cfg.sh> may be used to generate per sequence and per
rate point configuration files for a variety of common test conditions.

Parameters are set by the last value encountered on the command line.
Therefore if a setting is set via a configuration file, and then a
subsequent command line parameter changes that same setting, the command
line parameter value will be used.
