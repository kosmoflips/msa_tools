package GenKode::Commons;
use strict;
use warnings;

use Data::Dumper;
# shared small subs

use base 'Exporter';
our @EXPORT = qw/
chkfile
get_all_input_files
getorder
/;

our @EXPORT_OK = qw/

/;

our %EXPORT_TAGS = (
# translation=>[qw//]
);

sub get_all_input_files {
	my ($paths,$fromformats)=@_; # A ref, A ref
	my $filelist;
	return $filelist if !$paths;
	foreach my $f (@$paths) {
		if (-d $f) {
			my $dirf=getfile_dir($f,$fromformats);
			push @$filelist, @$dirf if $dirf;
		}
		elsif (chkfile($f,$fromformats)) {
			push @$filelist, $f;
		}
	}
	return $filelist;
}
sub getfile_dir { #get files from dir
	my ($dir,$fromformats)=@_;
	my $flist;
	opendir (my $dh, $dir);
	chdir $dir;
	while (my $f=readdir $dh) {
		next if $f=~/^\.+$/;
		if (chkfile($f,$fromformats)) {
			push @$flist, (File::Spec->catfile($dir,$f));
		}
	}
	return $flist;
}
sub chkfile {
	my ($file,$fromformats)=@_; #$fromformats is A ref
	my $code=0;
	if (-e $file and !-d $file and !-z $file and -T $file) {
		if ($fromformats and ref $fromformats=~/ARRAY/) {
			foreach my $fmt (@$fromformats) {
				if ($file=~/$fmt  $/xi) {
					$code=1;
					last;
				}
			}
		} else {
			$code=1;
		}
	}
	return $code;
}

sub getorder { # for reorder seq names. input a file of listed names, return A ref
	my ($sortfile)=@_;
	open (my $fh, $sortfile) or return undef;
	my $order;
	while (<$fh>) {
		next if !/\S/;
		next if /^#/;
		chomp;
		my (@c)=split /\t/; #in case this is a file for rename
		push @$order, $c[0];
	}
	return $order;
}

1;