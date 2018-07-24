
Format and processing of cfg/*.yaml by `gen-cfg.pl`
===================================================

All YAML config spec files are merged together (see merging rules) prior to
processing.  This has an important side-effect that two configuration
categories with the same name will be merged together rather than being
evaluated seperately.

The YAML config spec contains two top-level structures:

  `categories`, a set of configuration categories, each containing
  a set of sequence names and common encoder and decoder options.
  All sequences in a category will use the same common options.

  `sequences`, a set of sequence names, each describing common values
  used in the generation of the configuration files in all categories.
  The main use is to specify properties of the source data, such as its
  location file location, or processing options specific only to the
  sequence.

Configuration generation proceeds as follows:

 - For each category, a set of sequences are determined
 - For each sequence, a set of variants are determined
 - For each variant, configuration files are generated and written as:
    - `$prefix/$category/$sequence/$variant/encoder.cfg`
    - `$prefix/$category/$sequence/$variant/decoder.cfg`
    - `$prefix/$category/$sequence/$variant/pcerror.cfg`

To generate an encoder, decoder or pcerror configuration, cfgOptions are
gathered in the following order:

 - `encflags`/`decflags`/`pcerrorflags` from the global 'sequences'
    for the current sequence name
 - `encflags`/`decflags`/`pcerrorflags` for the current variant
    from the current category
 - `encflags`/`decflags`/`pcerrorflags` for the current variant
    from the current category-sequence
 - `encflags`/`decflags`/`pcerrorflags` set globally

Semantics of a yaml-cfg-file
----------------------------

The following description of the YAML config spec uses the following
conventions:

 - `.a=b` represents a map (associative array, dictionary, etc.,) with
    the key `a` and value `b`.  YAML in-line style: `{ a: b }`

 - `[]` represents a list of values. YAML in-line style: `[ ... ]`

 - `$x` represents a value (the value may also be a key in a map).

 - `/` represents the top-level YAML document.


### Definition of $cfgOption

A $cfgOption represents one of the following structures to generate a
configuration option in the form `$key: $value`

  - `.$key=$value` — General case
  - `.$key=[].$variant=.$value` — Applies only to the given $variant


### Top-level definitions

- `/.categories=.$categoryName=`...  
  A configuration category

- `/.sequences=.$sequenceName=`...  
  A global sequence definition

- `/.sequence-base-dir=$value`  
  A global base directory that may be overriden by a sequences' 
  `.base-dir` and `.base-norm-dir` values.

- `/.pcerrorflags=[].$cfgOption`  
  (optional) an ordered list of global options for pcerror.cfg

- `/.encflags=[].$cfgOption`  
  (optional) an ordered list of global options for encoder.cfg

- `/.decflags=[].$cfgOption`  
  (optional) an ordered list of global options for decoder.cfg

### Inside `/.sequences=.$sequenceName=`...

- `.src=$value`  
  The source PLY filename for encoding

- `.src-dir=$value`  
  (optional) The directory name containing the .src file

- `.base-dir=$value`  
  (optional) A path to a directory containing .src-dir

- `.norm=$value`  
  (optional) The source PLY filename with normals data

- `.norm-dir=$value`  
  (optional) The directory containing the .norm file

- `.base-norm-dir=$value`  
  (optional) A path to a directory containint `.norm-dir`

- `.pcerrorflags=[].$cfgOption`  
  (optional) an ordered list of sequence-global options for pcerror.cfg

- `.encflags=[].$cfgOption`  
  (optional) an ordered list of sequence-global options for encoder.cfg

- `.decflags=[].$cfgOption`  
  (optional) an ordered list of sequence-global options for decoder.cfg

### Inside `/.categories=.$categoryName=`...

- `.encflags=[].$cfgOption`  
  (optional) an ordered list of category-specific options for encoder.cfg

- `.decflags=[].$cfgOption`  
  (optional) an ordered list of category-specific options for decoder.cfg

- `.pcerrorflags=[].$cfgOption`  
  (optional) an ordered list of category-specific options for pcerror.cfg

- `.sequences=...`  
  A set of sequences to generate configurations for in the context
  of the current category

### Inside `/.categories=.$categoryName=.sequences=.$sequenceName=...`

- `.encflags=[].$cfgOption`  
  (optional) an ordered list of category-sequence-specific options
  for encoder.cfg

- `.decflags=[].$cfgOption`  
  (optional) an ordered list of category-sequence-specific options
  for decoder.cfg

- `.pcerrorflags=[].$cfgOption`  
  (optional) an ordered list of category-sequence-specific options
  for pcerror.cfg

## Merging rules

Multiple YAML config spec files are recursively merged as follows:

- src:*      → dst:undef  ⇒ assign src to dst
- src:scalar → dst:scalar ⇒ replaced
- src:hash   → dst:hash   ⇒ recursive merge of key-value pairs
- src:list   → dst:scalar ⇒ assign [src, dst] to dst
- src:list   → dst:list   ⇒ assign [src, dst] to dst
