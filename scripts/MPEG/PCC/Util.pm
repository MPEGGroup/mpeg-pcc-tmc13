package MPEG::PCC::Util;
use strict;
use warnings;

use Exporter qw(import);
our @EXPORT_OK = qw(uniq);

##
# From List::MoreUtils::uniq, licensed as follows:
# > This library is free software; you can redistribute it and/or modify
# > it under the same terms as Perl itself, either Perl version 5.8.4
# > or, at your option, any later version of Perl 5 you may have
# > available.
sub uniq (@)
{
	my %seen = ();
	my $k;
	my $seen_undef;
	grep { defined $_ ? not $seen{$k = $_}++ : not $seen_undef++ } @_;
}

1;
