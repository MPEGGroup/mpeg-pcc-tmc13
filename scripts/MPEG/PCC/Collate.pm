package MPEG::PCC::Collate;
use utf8;
use strict;
use warnings;
use MPEG::PCC::Util qw{uniq};
use List::Util qw{min max pairmap pairgrep};
use Scalar::Util qw{looks_like_number};

use Exporter qw(import);
our @EXPORT = qw(
	accumulateResult summariseResult
);

our @out_order_cols = qw{
	config sequence variant framenum
	enc.status dec.status dec.enc.match
	enc.wtime enc.utime
	dec.wtime dec.utime
	enc.ext.utime dec.ext.utime enc.ext.maxrssk dec.ext.maxrssk

	enc.ext.bits enc.bits enc.ext.bpp enc.bpp
	enc.bpp.geometry enc.bpp.colour enc.bpp.reflectance
	enc.bits.geometry enc.bits.colour enc.bits.reflectance
	enc.bits.frameindexs

	src.numpoints dec.numpoints
	src.framecount dec.framecount

	dec.y-psnr dec.cb-psnr dec.cr-psnr dec.reflectance-psnr
	dec.post-recolour.y-psnr
	dec.post-recolour.cb-psnr dec.post-recolour.cr-psnr
	dec.post-recolour.reflectance-psnr
	dec.y-hpsnr dec.cb-hpsnr dec.cr-hpsnr dec.reflectance-hpsnr
	dec.d1-psnr dec.d2-psnr
};

our @fold_mean_geom = qw{
};

our @fold_mean_arith = qw{
	enc.bpp
	enc.bpp.geometry enc.bpp.colour enc.bpp.reflectance
	dec.y-psnr dec.cb-psnr dec.cr-psnr dec.reflectance-psnr
	dec.d1-psnr dec.d2-psnr
	dec.post-recolour.y-psnr
	dec.post-recolour.cb-psnr dec.post-recolour.cr-psnr
	dec.post-recolour.reflectance-psnr
};

our @fold_sum = (@fold_mean_arith, qw{
	src.framecount dec.framecount
	enc.ext.utime enc.wtime enc.utime
	dec.ext.utime dec.wtime dec.utime
	enc.ext.bits
	enc.bits enc.bits.geometry enc.bits.colour enc.bits.reflectance
	enc.bits.frameindexs
	src.numpoints dec.numpoints
});

our @fold_max = qw{
	enc.ext.maxrssk dec.ext.maxrssk
	dec.y-hpsnr dec.cb-hpsnr dec.cr-hpsnr dec.reflectance-hpsnr
};

our @fold_uniq = qw{ enc.status dec.status dec.enc.match };

##
# accumulate results for frame into %$result_acc
#
sub accumulateResult {
	my ($result_acc, $frame) = @_;

	# accumulate additively
	foreach my $key (@fold_sum) {
		next unless exists $frame->{$key} && $frame->{$key} ne '';
		$$result_acc{"count.$key"} = 1 + ($$result_acc{"count.$key"} || 0);
		$$result_acc{$key} = $frame->{$key} + ($$result_acc{$key} || 0);
		delete $frame->{$key};
	}

	# accumulate maximum
	foreach my $key (@fold_max) {
		next unless exists $frame->{$key};
		next unless $frame->{$key};
		$$result_acc{$key} = max($frame->{$key},$$result_acc{$key} || 0);
		delete $frame->{$key};
	}

	# accumulate log for geometric mean calculation
	foreach my $key (@fold_mean_geom) {
		next unless exists $frame->{$key};
		if (!$frame->{$key} > 0) {
			$$result_acc{"invalid.$key"} = 1;
			next;
		}
		$$result_acc{"count.$key"} = 1 + ($$result_acc{"count.$key"} || 0);
		$$result_acc{$key} = log($frame->{$key}) + ($$result_acc{$key} || 0);
		delete $frame->{$key};
	}

	# accumulate unique values
	foreach my $key (@fold_uniq) {
		next unless exists $frame->{$key};
		$$result_acc{$key} = [uniq($frame->{$key}, @{$$result_acc{$key}})];
		delete $frame->{$key};
	}

	$result_acc;
}

##
# reduce the accumulation state to summary results
#
sub summariseResult {
	my ($result_acc) = @_;

	# geometric means (relies on sum of logs)
	foreach my $key (@fold_mean_geom) {
		next unless exists $$result_acc{$key};
		if ($$result_acc{"invalid.$key"}) {
			$$result_acc{$key} = 'NaN';
			next;
		}
		$$result_acc{$key} =
			 exp($$result_acc{$key} / $$result_acc{"count.$key"});
	}

	# arithmetic mean
	foreach my $key (@fold_mean_arith) {
		$$result_acc{$key} = !$$result_acc{"count.$key"} ? undef :
			$$result_acc{$key} / $$result_acc{"count.$key"};
	}

	# convert to friendly string
	foreach my $key (@fold_uniq) {
		$$result_acc{$key} = join ':', @{$$result_acc{$key}}
			if defined $$result_acc{$key};
	}

	# ratio from totals
	$$result_acc{'enc.ext.bpp'} =
		$$result_acc{'enc.ext.bits'} / $$result_acc{'src.numpoints'}
		if $$result_acc{'src.numpoints'};

	# tidyup any formatting
	map {
		if (looks_like_number($_)) {
			$_ = sprintf "%f", $_;
			s/\.?0*$//;
		}
	} values %$result_acc;

	return $result_acc;
}

1;
