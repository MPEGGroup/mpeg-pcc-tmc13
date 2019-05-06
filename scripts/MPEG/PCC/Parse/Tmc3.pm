package MPEG::PCC::Parse::Tmc3;
use strict;
use warnings;

use Exporter qw(import);
our @EXPORT = qw(readEncLog readDecLog);

##
# parse output of encoder log
sub readEncLog {
	my ($file) = @_;

	open my $fh, '<', $file
		or return {};

	my %result;
	while (<$fh>) {
		chomp;
		if (m{^(positions|colors|reflectances|\w+) bitstream size (\d+) B \((\d+(\.\d+(e[+-]\d+)?)?) bpp\)}) {
			my %map = (
				positions => 'geometry',
				colors => 'colour',
				reflectances => 'reflectance',
			);
			my $key = $map{$1} || $1;

			no warnings;
			$result{"enc.bits.$key"} += $2 * 8;
			$result{"enc.bpp.$key"} += $3;
			next;
		}

		if (m{^Total bitstream size (\d+) B}) {
			$result{'enc.bits'} = $1 * 8;
			next;
		}

		if (m{^Processing time \(wall\): (\d+(\.\d+)?) s}) {
			$result{'enc.wtime'} = $1;
			next;
		}

		if (m{^Processing time \(user\): (\d+(\.\d+)?) s}) {
			$result{'enc.utime'} = $1;
			next;
		}
	}

	return \%result;
}

##
# parse output of decoder log
sub readDecLog {
	my ($file) = @_;

	open my $fh, '<', $file
		or return {};

	my %result;
	while (<$fh>) {
		chomp;
		if (m{^Processing time \(wall\): (\d+(\.\d+)?) s}) {
			$result{'dec.wtime'} = $1;
			next;
		}

		if (m{^Processing time \(user\): (\d+(\.\d+)?) s}) {
			$result{'dec.utime'} = $1;
			next;
		}
	}

	return \%result;
}

1;
