package MPEG::PCC::CSV;

use strict;

use Text::CSV;

use Exporter qw(import);
our @EXPORT = qw(
	LoadFile
);

sub LoadFile {
	my ($name) = @_;
	my $fh;

	if (ref $name eq 'GLOB') {
		$fh = $name;
	} else {
		open $fh, "<:encoding(utf8)", $name or die "$name: $!";
	}

	# drop any comment lines at the start
	my $comment = "";
	my $first_non_comment_line;
	while (my $line = <$fh>) {
		unless ($line =~ m{^\#}) {
			# Avoiding dependency from IO::Unread::unread $fh, $line;
			$first_non_comment_line = $line;
			last;
		}
		$comment .= $line;
	}

	my $csv = Text::CSV->new({binary => 1});
	my $cols = do {
		open my $tmpfh, '<', \$first_non_comment_line or die;
		$csv->getline($tmpfh);
	};

	# handle empty file case
	return wantarray ? ([], [], $comment) : [] unless defined $cols;

	for (my $i = 0; $i < scalar @$cols; $i++) {
		# add fake column names if any are missing
		next if defined $cols->[$i] && !($cols->[$i] =~ m{^\s*$});
		$cols->[$i] = "_$i";
	}
	$csv->column_names($cols);
	my $rows = $csv->getline_hr_all($fh);

	return wantarray ? ($rows, $cols, $comment) : $rows;
}

1;
