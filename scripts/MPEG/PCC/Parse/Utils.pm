package MPEG::PCC::Parse::Utils;
use strict;
use warnings;

use Exporter qw(import);
our @EXPORT = qw(readFileFirstLine);

##
# cat the first line of a file.
sub readFileFirstLine {
	my ($file) = @_;

	open my $fh, '<', $file
		or return ();
	chomp (my $line = <$fh>);
	return $line;
}

1;
