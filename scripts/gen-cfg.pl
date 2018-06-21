#!/usr/bin/perl

use Digest::MD5;
use File::Path qw(make_path);
use File::Basename qw(basename);
use Getopt::Long;
use List::MoreUtils;
use Module::Load::Conditional qw(can_load);
use Pod::Usage;
use YAML;
use strict;

=head1 NAME

gen-cfg.pl - Generate experiment configuration from yaml specification

=head1 SYNOPSIS

gen-cfg.pl [options] [yaml-config-spec ...]

=head1 OPTIONS

=over 4

=item B<--prefix>=dir

Sets the output path for the generated configuration tree.

=item B<--output-src-glob-sh>

=item B<--no-output-src-glob-sh> (default)

Do not generate src-glob.sh files describing source locations

=item B<--skip-sequences-without-src> (default)

=item B<--no-skip-sequences-without-src>

Do not generate configuration files for sequences that have an empty 'src'
field in the yaml specification.  This option is permits a later yaml spec
to effectively remove a sequence from being used in an experiment.

It may be useful to disable this option when generating config files when
the source location of the input data is not known.

=item B<--batch-template>=script.tt

Generate a script from the template script.tt.  The output is written to
the current working directory.

=item B<--experiment-name>=name

A generic name that may be referenced by a batch template for any purpose.

=back

=head1 Config specification files



=cut


##
# Command line processing
my $do_help = '';
my $output_src_glob_sh = 0;
my $skip_sequences_without_src = 1;
my $experiment_name = '';
my $batch_template = '';
my $prefix = '.';
GetOptions(
	'help' => \$do_help,
	'prefix=s' => \$prefix,
	'output-src-glob-sh!' => \$output_src_glob_sh,
	'skip-sequences-without-src!' => \$skip_sequences_without_src,
	'batch-template=s' => \$batch_template,
	'experiment-name=s' => \$experiment_name,
);

##
# display help text and exit if asked, or if no config is provided
pod2usage(0) if $do_help;
pod2usage(1) unless @ARGV;

##
# load all yaml snippets and merge into a single description
#
my @origins = @ARGV;
my %cfg;
while (@ARGV) {
	my $fname = shift @ARGV;
	my $cfg = YAML::LoadFile($fname) or die "$fname: $!";
	merge(\%cfg, $cfg);
}

##
# dump the merged configuration (allows reproduction)
#
YAML::DumpFile("$prefix/config-merged.yaml", \%cfg);

##
# generate encoder/decoder configuration files
#

# list of configured jobs
my @jobs;

# this just makes later code look simpler
my $cfg = \%cfg;

# iterate over each configuration and described sequences
foreach my $cat_name (sort keys %{$cfg->{categories}}) {
	my $cat = $cfg->{categories}{$cat_name};

	foreach my $seq_name (sort keys %{$cat->{sequences}}) {
		my $cat_seq = $cat->{sequences}{$seq_name};
		my $seq = $cfg->{sequences}{$seq_name};

		unless (exists $seq->{gops}) {
			genSeqVariants($cat, $cat_name, $cat_seq, $seq_name, $seq, $seq);
			next;
		}

		# split sequence into groups of pictures for parallel execution
		my $gop_idx = 0;
		foreach my $gop (@{$seq->{gops}}) {
			my $gop_idx = sprintf "%03d", $gop_idx++;
			my $gop_name = "${seq_name}_gop${gop_idx}";
			genSeqVariants($cat, $cat_name, $cat_seq, $gop_name, $gop, $seq);
		}
	}
}

##
# write out batch-job configuration
#
if ($batch_template) {
	can_load(modules => {'Template' => undef})
		|| die $Module::Load::Conditional::ERROR;

	my $output = basename($batch_template,'.tt');

	my $vars = {
		jobs => \@jobs,
		experiment_name => $experiment_name,
	};

	my $tt = Template->new({
		RELATIVE => 1, # allow relative paths as $batch_template
		ABSOLUTE => 1, # allow absolute paths too
	}) || die "$Template::ERROR\n";
	$tt->process($batch_template, $vars, $output) || die $tt->error(), "\n";
}


sub genSeqVariants {
	my ($cat, $cat_name, $cat_seq, $seq_name, $gop, $seq) = @_;

	# if sequence source isn't defined at top level, skip
	if ($skip_sequences_without_src) {
		next unless defined $gop->{src};
	}

	# generate the list of variants (if any)
	my @variants = List::MoreUtils::uniq (
		# $cat.sequences.$name.$variant:
		(grep {
				my $ref = $cat_seq->{$_};
				ref $ref eq 'HASH'
				and (exists $ref->{encflags} || exists $ref->{decflags})
			} keys %$cat_seq),

		# $cat.sequences.$name.encflags[].$param.$variant:
		variants_from_node($cat_seq->{encflags}),

		# $cat.encflags[].$param.$variant
		variants_from_node($cat->{encflags}),

		# $seq.$name.encflags[].$param.$variant
		variants_from_node($seq->{encflags}),
	);

	# handle the case of no variants: single case with defaults
	push @variants, undef unless @variants;

	# for each variant, derive the encoder options
	#   NB: in the case of no variants, $var = undef
	foreach my $var (sort @variants) {

		my $cfgdir =
			join '/', grep {defined} ($prefix,$cat_name,$seq_name,$var);
		print "$cfgdir\n";
		make_path($cfgdir);
		push @jobs, "$cfgdir/";

		##
		# input sequence file name
		if ($gop->{src} && $output_src_glob_sh) {
			my $src_seq = join '/', grep {defined} (
				(List::MoreUtils::firstval {defined}
					$seq->{'base-dir'},
					$cfg->{'sequence-base-dir'}),
				$seq->{'src-dir'},
				$gop->{src},
			);

			open my $fd, ">", "$cfgdir/src-glob.sh";
			print $fd "$src_seq\n";
		}

		if ($gop->{norm} && $output_src_glob_sh) {
			my $norm_seq = join '/', grep {defined} (
				(List::MoreUtils::firstval {defined}
					$seq->{'base-norm-dir'},
					$seq->{'base-dir'},
					$cfg->{'sequence-base-dir'}),
				$seq->{'norm-dir'},
				$gop->{norm},
			);

			open my $fd, ">", "$cfgdir/norm-glob.sh";
			print $fd "$norm_seq\n";
		}

		##
		# encoder configuration
		my @encflags = (
			params_from_node($seq->{encflags}),
			params_from_node($cat->{encflags}, $var),
			params_from_node($cat_seq->{encflags}, $var),
			params_from_node($cat_seq->{$var}{encflags}),
		);

		# evaluate any value expressions
		eval_exprs(\@encflags, $cat_seq, $seq);

		write_cfg("$cfgdir/encoder.cfg", \@encflags);

		##
		# decoder configuration
		my @decflags = (
			params_from_node($seq->{decflags}),
			params_from_node($cat->{decflags}, $var),
			params_from_node($cat_seq->{decflags}, $var),
			params_from_node($cat_seq->{$var}{decflags}),
		);

		# evaluate any value expressions
		eval_exprs(\@decflags, $cat_seq, $seq);

		write_cfg("$cfgdir/decoder.cfg", \@decflags);

		##
		# pcerror configuration
		my @pcerrorflags = (
			params_from_node($seq->{pcerrorflags}),
			params_from_node($cat->{pcerrorflags}),
			params_from_node($cat_seq->{pcerrorflags}, $var),
			params_from_node($cat_seq->{$var}{pcerrorflags}),
		);

		# evaluate any value expressions
		eval_exprs(\@pcerrorflags, $cat_seq, $seq);

		write_cfg("$cfgdir/pcerror.cfg", \@pcerrorflags) if (@pcerrorflags);
	}
}


#############################################################################
# utilities

##
# keywise merge $src into $dst, following the following merge rules:
#  - *      -> undef  = copy
#  - scalar -> scalar = replace
#  - hash   -> hash   = recurse
#  - list   -> scalar = merge unique items (scalars only)
#  - list   -> list   = merge unique items (scalars only)
sub merge {
	my ($dst, $src) = @_;

	unless (defined $dst) {
		$$dst = $$src;
		return;
	}

	##
	# overwrite existing scalar
	unless (ref $src) {
		$$dst = $src;
		return;
	}

	if (ref $src eq 'HASH') {
		foreach my $key (keys %$src) {
			##
			# copy sub-tree if key does not exist
			unless (exists $$dst{$key}) {
				$$dst{$key} = $$src{$key};
				next;
			}

			##
			# recurse to merge sub-tree
			if (ref $$dst{$key}) {
				merge($$dst{$key}, $$src{$key});
			}
			else {
				merge(\$$dst{$key}, $$src{$key});
			}
		}
		return;
	}

	##
	# merge arrays
	#  -- this is really only for an array of scalars
	if (ref $src eq 'ARRAY') {
		my @vals;
		push @vals, $$dst if ref $dst eq 'SCALAR';
		push @vals, @$dst if ref $dst eq 'ARRAY';
		push @vals, @$src;

		$$dst = [List::MoreUtils::uniq(@vals)] if ref $dst eq 'SCALAR';
		@$dst = List::MoreUtils::uniq(@vals) if ref $dst eq 'ARRAY';
		return;
	}
}


sub variants_from_node {
	my ($node) = @_ or return ();
	return
		map {keys %$_}
		grep {ref $_ eq 'HASH'}
		map {values %$_}
		grep {ref $_ eq 'HASH'}
		map {ref $_ eq 'ARRAY' ? @$_ : $_}
		@{$node};
}

sub params_from_node {
	my ($node, $variant) = @_;
	return () unless $node;

	my @params;
	my @todo = @$node;
	while (my $item = shift @todo) {
		# an unformatted string (not key:value)
		unless (ref $item) {
			push @params, [$item];
			next;
		}
		if (ref $item eq 'HASH') {
			while (my ($key, $value) = each %$item) {
				unless (ref $value) {
					# key:value without variants
					push @params, [$key, $value];
					next;
				}
				if (ref $value eq 'HASH') {
					# key:value with variants
					push @params, [$key, $value->{$variant}]
						if exists $value->{$variant};
					next;
				}
				warn "unhandled node for $value";
			}
		}
		if (ref $item eq 'ARRAY') {
			unshift @todo, @$item;
			next;
		}
	}
	return @params;
}

##
# Expand in-place all variables and eval statements in the list @$params,
# searching for substitutions as members of @context
sub eval_exprs {
	my ($params, @context) = @_;

	map {
		$_->[1] = eval_expr($_->[1], \@context) if exists $_->[1];
	} @$params;
}

##
# Return the exansion of $str given the context of the maps in @$contexts.
# Any variable subsitutions found in $str are searched in order
# of @$contexts.
sub eval_expr {
	my ($str, $contexts) = @_;

	# first find all variables and substitute their values
	while ($str =~ m/\$\{([^}]+)\}/gc) {
		my $var = $1;
		my $var_start = $-[0];
		my $var_len  = $+[0] - $-[0];
		foreach my $ctx (@$contexts) {
			next unless exists $ctx->{$var};
			substr $str, $var_start, $var_len, $ctx->{$var};
			pos $str = $var_start + length($ctx->{$var});
			last;
		}
	}

	# finally evaluate any eval expressions
	pos $str = 0;
	while ($str =~ m/\$eval\{([^}]+)\}/gc) {
		my $expr = $1;
		my $expr_start = $-[0];
		my $expr_len  = $+[0] - $-[0];
		my $val = eval "$expr";
		substr $str, $expr_start, $expr_len, $val;
		pos $str = $expr_start + length($val);
	}

	return $str;
}


##
# Print configuration @$opts, to $fd; with one entry per line and where
# each entry in @$opts is either a [key, value] pair to be joined with
# ": ", or just [key].
sub print_cfg {
	my ($fd, $opts) = @_;

	print $fd "# This file was automatically generated from:\n";
	print $fd "#   $_\n" foreach (@origins);

	local $\ = "\n";
	foreach my $opt (@$opts) {
		print $fd join(": ", @$opt);
	}
}


##
# print config to file iff it differs from file's contents.
# (ie, don't touch mtime if unchanged)
sub write_cfg {
	my ($filename, $flags) = @_;

	# format config in memory
	my $new_cfg = "";
	open my $fd, ">:encoding(utf8)", \$new_cfg;
	print_cfg($fd, $flags);
	close $fd;

	# hash it
	my $md5_new = Digest::MD5->new;
	$md5_new->add($new_cfg);

	my $md5_old = Digest::MD5->new;
	if (-f $filename) {
		open $fd, "<:encoding(utf8)", $filename;
		$md5_old->addfile($fd);
		close $fd;
	}

	if ($md5_new->digest ne $md5_old->digest) {
		print "writing $filename\n";
		open $fd, ">:bytes", $filename;
		print $fd $new_cfg;
		close $fd;
	}
}
