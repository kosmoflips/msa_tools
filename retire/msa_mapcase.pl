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

#input 2 msa files, one aligned but no case difference. one may not be aligned but have correct cases


my ($help,$infile,$reffile);
GetOptions(
	'help'=>\$help,
	'infile|file=s{1}'=>\$infile,
	'reffile=s{1}'=>\$reffile
);


if (!$infile or !$reffile or !-e $infile or !-e $reffile) {
	print <<HELP;
--------------------------------------
[-i INPUT_FILE] #alignment ok but lacks correct case
[-r REF_FILE] #has correct case
# sequences in both files must have the same name (case insensitive) in order to be case-remapped
--------------------------------------
HELP
exit;
}

my $mref=GenKode::MSA->new;
my $m=GenKode::MSA->new;
$mref->import($reffile);
$m->import($infile);

my $refname;
foreach my $id (1..$mref->lastid) {
	$refname->{lc $mref->getname($id)}=$id;
}

foreach my $id (1..$m->lastid) {
	if ($refname->{lc $m->getname($id)}) {
		my $rseq=$mref->getseq($refname->{lc $m->getname($id)});
		$m->map_case($id,$rseq);
	}
}

my $ofile=$infile.'_casemapped.fst';
$m->export($ofile, 'fst',1);

printf "- output write to %s\n", $ofile;