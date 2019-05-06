package MPEG::PCC::Parse::PcError;
use strict;
use warnings;

use Exporter qw(import);
our @EXPORT = qw(readDistortion);

##
# mapping table for readDistortion
our %readDistortion_key2key = (
	'h.        (p2point)' => 'd1-hmse',  # hausdorff error
	'h.,PSNR   (p2point)' => 'd1-hpsnr', # hausdorff error
	'h.        (p2plane)' => 'd2-hmse',  # hausdorff error
	'h.,PSNR   (p2plane)' => 'd2-hpsnr', # hausdorff error
	'mseF      (p2point)' => 'd1-mse',
	'mseF,PSNR (p2point)' => 'd1-psnr',
	'mseF      (p2plane)' => 'd2-mse',
	'mseF,PSNR (p2plane)' => 'd2-psnr',
	'c[0],    F'          => 'y-mse',
	'c[1],    F'          => 'cb-mse',
	'c[2],    F'          => 'cr-mse',
	'c[0],PSNRF'          => 'y-psnr',
	'c[1],PSNRF'          => 'cb-psnr',
	'c[2],PSNRF'          => 'cr-psnr',
	'r,       F'          => 'reflectance-mse',
	'r,PSNR   F'          => 'reflectance-psnr',
	'h.c[0],    F'        => 'y-hmse',
	'h.c[1],    F'        => 'cb-hmse',
	'h.c[2],    F'        => 'cr-hmse',
	'h.c[0],PSNRF'        => 'y-hpsnr',
	'h.c[1],PSNRF'        => 'cb-hpsnr',
	'h.c[2],PSNRF'        => 'cr-hpsnr',
	'h.r,       F'        => 'reflectance-hmse',
	'h.r,PSNR   F'        => 'reflectance-hpsnr',
);

##
# parse output of pc_error
sub readDistortion {
	my ($file, $key_prefix) = @_;

	open my $fh, '<', $file
		or return {};

	my %result;

	# skip over all the preamble,
	while (<$fh>) {
		if (m{^PCC quality measurement software, version (.*)}) {
			$result{"$key_prefix.dmetric.version"} = $1;
			next;
		}
		last if m{^3. Final \(symmetric\).};
	}

	# read in the record
	our %readDistortion_key2key;
	while (<$fh>) {
		chomp;
		# change in indentation breaks the record
		last unless m{^ };
		s/^\s*//;
		my ($key, $val) = split /\s*:\s*/;
		$result{$key_prefix.$readDistortion_key2key{$key}} = $val;
	}

	return \%result;
}

1;
