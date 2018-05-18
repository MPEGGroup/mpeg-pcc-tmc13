
# Warning

The following sub-directories containing configuration files are
automatically generated from the .yaml master files and should not
be manually edited:

 lossless-geom-lossy-reflectance
 lossless-geom-lossy-texture
 lossless-geom-no-attrs
 lossy-geom-lossy-attrs
 lossy-geom-no-attrs

Any manual changes that are not reflected in the .yaml master files
is likely to be lost or to cause re-integration issues.

# How to generate the per-data point config files

Run the `../scripts/gen-cfg.sh` from within the cfg directory.
