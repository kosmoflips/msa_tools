#!/usr/bin/perl
use strict;
use warnings;

use File::Spec;
use File::Path;
use File::Copy;
use File::Temp;
use Storable qw/:DEFAULT nstore dclone/;
use Data::Dumper;
use Getopt::Long;

use lib '.';
use GenKode::MSA;
use GenKode::Commons;

# join 2 msa files . must name match


my ($help,@files,$addblock,$ref);
GetOptions(
	'help'=>\$help,
	'add'=>\$addblock,
	'ref=i'=>\$ref,
	'file=s{2,}'=>\@files,
);

my $filelist=GenKode::Commons::get_all_input_files(\@files);

if (!$filelist or $help) {
	print <<HELP;
--------------------------------------
[-f FILE1 2 ...] #at least 2 files. file paths in order from 5' to 3'
[-a] #optional. separate two blocks using a string of "N"
[-r] #optional. specify the Xth file given is the base file
--------------------------------------
HELP
exit;
}

if (!$ref or $ref>scalar @$filelist) {
	$ref=1;
}

my $m1=GenKode::MSA->new;
$m1->import($filelist->[$ref-1]);


my $c=0;
foreach my $file (@$filelist) {
	$c++;
	printf "- %s", $file;
	if ($c==$ref) {
		print " (root file)\n";
		next;
	} else {
		printf "\n";
	}
	my $m2=GenKode::MSA->new;
	my $readchk=$m2->import($file);
	if (!$readchk) {
		print "  can't parse file. skip.\n";
		next;
	}
	my $refname;
	foreach my $id (1..$m2->lastid) {
		$refname->{lc $m2->getname($id)}=$id;
	}

	foreach my $id (1..$m1->lastid) {
		my $name1=$m1->getname($id);
		if ($refname->{lc $name1}) { #has name match, do conj
			my $lseq=$m1->getseq($id);
			my $lseqlen=length $lseq;
			my $rseq=$m2->getseq($refname->{lc $name1});
			my $rseqlen=length $lseq;
			my $sc=$lseq;
			if ($lseqlen<$m1->maxlen) {
				$lseq.='.' x ($m1->maxlen - $lseqlen);
			}
			if ($rseqlen<$m2->maxlen) {
				$rseq.='.' x ($m2->maxlen - $rseqlen);
			}
			if ($ref>$c) { #put before
				$sc=$rseq.($addblock?'n' x 15:'').$sc;
			} else { #put after
				$sc.=($addblock?'n' x 15:'').$rseq;
			}
			$m1->{$id}{seq}=$sc;
		} else { #fill gap
			my $sc0=$m1->{$id}{seq};
			my $add='.' x ($m2->maxlen);
			if ($ref>$c) { #put before
				$sc0=$add.($addblock?'n' x 15:'').$sc0;
			} else { #put after
				$sc0.=($addblock?'n' x 15:'').$add;
			}
			$m1->{$id}{seq}=$sc0;
		}
	}
}

my $ofile=$filelist->[$ref-1].'_joined.fst';
$m1->export($ofile, 'fst',1);

printf "\n-- output write to %s\n", $ofile;