package MPEG::PCC::Parse::Time;
use strict;
use warnings;

use Exporter qw(import);
our @EXPORT = qw(readTime);

##
# parse output of /bin/time.
# returns (user_time, maxrss)
sub readTime {
	my ($file) = @_;

	open my $fh, '<', $file
		or return ('?');
	chomp (my $line = <$fh>);

	my $utime;
	my $maxrssk;
	foreach (split / /, $line) {
		if (m{^(\d+\.\d+)user$})  { $utime = $1; next; }
		if (m{^(\d+)maxresident}) { $maxrssk = $1; next; }
	}

	return ($utime, $maxrssk);
}

1;
