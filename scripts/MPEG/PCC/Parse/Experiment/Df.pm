package MPEG::PCC::Parse::Experiment::Df;
use strict;
use warnings;

use MPEG::PCC::Parse::Tmc3;
use MPEG::PCC::Parse::Time;
use MPEG::PCC::Parse::Utils;
use MPEG::PCC::Parse::Ply;
use MPEG::PCC::Parse::PcError;

use Exporter qw(import);
our @EXPORT = qw(
	readTmc3Results readTmc3ResultsOneFrame readTmc3ResultsOneBin
);

##
# One frame per file mode (ie not sequence encoding)
sub readTmc3Results {
	my ($base_path, $src_file) = @_;
	my ($frame) = $src_file =~ m{/([^/]+)$};

	my $ret_ply = readTmc3ResultsOneFrame(@_);
	my $ret_bin = readTmc3ResultsOneBin(@_);

	return {%$ret_ply, %$ret_bin, frame => $frame};
}


##
# process a binary, allowing results to be aggregated
sub readTmc3ResultsOneBin {
	my ($base_path) = @_;

	my $file_bytes = -s "$base_path.bin";
	my ($enc_utime, $enc_maxrssk) = readTime("$base_path.bin.time");
	my ($dec_utime, $dec_maxrssk) = readTime("$base_path.bin.decoded.time");

	my $enc_log = readEncLog("$base_path.bin.log");
	my $dec_log = readDecLog("$base_path.bin.decoded.log");

	my $enc_status = readFileFirstLine("$base_path.bin.status");
	my $dec_status = readFileFirstLine("$base_path.bin.decoded.status");

	my %ret = (
		"enc.ext.bits" => $file_bytes * 8,
		"enc.ext.utime" => $enc_utime,
		"dec.ext.utime" => $dec_utime,
		"enc.ext.maxrssk" => $enc_maxrssk,
		"dec.ext.maxrssk" => $dec_maxrssk,
		"enc.status" => $enc_status,
		"dec.status" => $dec_status,
		%$enc_log,
		%$dec_log,
	);

	return \%ret;
}

##
# process a frame, allowing results to be aggregated
sub readTmc3ResultsOneFrame {
	my ($base_path, $src_file) = @_;

	my ($num_src_points) = readPly("$src_file");
	my ($num_dec_points) = readPly("$base_path.bin.decoded.ply");

	my $distortion_e2e =
		readDistortion("$base_path.bin.decoded.pc_error", "dec.");

	my $distortion_prc =
		readDistortion(
			"$base_path.bin.decoded.pc_error_postrecolour",
			"dec.post-recolour.");

	my ($enc_md5,undef) =
		split / /, readFileFirstLine("$base_path.bin.ply.md5") // "";

	my ($dec_md5,undef) =
		split / /, readFileFirstLine("$base_path.bin.decoded.ply.md5") // "";

	my $dec_enc_match = "mismatch";
	if (!$enc_md5 || !$dec_md5) {
		$dec_enc_match = "missing";
	}
	elsif ($enc_md5 eq $dec_md5) {
		$dec_enc_match = "ok"
	}

	my %ret = (
		"src.numpoints" => $num_src_points,
		"src.framecount" => 1,
		"dec.numpoints" => $num_dec_points,
		"dec.framecount" => 1,
		"dec.enc.match" => $dec_enc_match,
		%$distortion_e2e,
		%$distortion_prc,
	);

	return \%ret;
}
