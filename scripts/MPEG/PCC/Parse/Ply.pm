package MPEG::PCC::Parse::Ply;
use strict;
use warnings;

use Exporter qw(import);
our @EXPORT = qw(readPly);

##
# cache of ply data to reduce number of lookups
our %ply_cache;

##
# parse ply file for interesting parameters
sub readPly {
	my ($file) = @_;

	our %ply_cache;
	return $ply_cache{$file} if exists $ply_cache{$file};

	open my $fh, '<', $file
		or return ();

	# check it is a ply file
	return () unless (<$fh> =~ m{^ply});

	my $numpoints;
	while (<$fh>) {
		chomp;
		# avoid the data section
		last if m{^end_header};

		if (m{^element vertex (\d+)}) { $numpoints = $1; next; }
	}

	return $ply_cache{$file} = $numpoints;
}

1;
