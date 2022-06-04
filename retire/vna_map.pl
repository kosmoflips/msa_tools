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

# input ALN, input foldings, add folding to ALN. by name match, case insensitive
use lib '.';
use GenKode::Fold;
use GenKode::Commons;

my ($help, $alnfile, @foldfile);
GetOptions(
	'help'=>\$help,
	'file|aln=s{1}'=>\$alnfile,
	'vna=s{1,}'=>\@foldfile, #can be multiple in case of many ct's
);

if ($help or !$alnfile or !@foldfile) {
	die <<HELP;

-------------------------------------
# add vna folding info into aligned sequences

>> input <<
[-f ALN_FILE] # contains aligned sequences
[-v Fold1.ct Fold2.vna Fold3_5.fst ...] # contain foldings for seq's in ALN_FILE.
  # match by name, case insensitive
  # if the same name appear multiple times, the LAST one will be actually mapped
-------------------------------------

HELP
}

my $foldfilelist=GenKode::Commons::get_all_input_files(\@foldfile);
my $v=GenKode::Fold->new;
foreach my $ff (@$foldfilelist) {
	$v->import($ff);
}
my $namemap;
foreach my $id (1..$v->lastid) {
	my $name=$v->getname($id);
	next if !$name;
	$namemap->{lc $name}=$id;
}

my $aln=GenKode::Fold->new;
$aln->import($alnfile);
for my $id (1..$aln->lastid) {
	my $name=$aln->getname($id);
	my $vid=$namemap->{lc $name}||0;
	if ($vid) {
		my $fold=GenKode::Fold::mkfold($v->{$vid}{seq},$v->{$vid}{vna});
		$aln->map_fold($id,$fold);
	}
}

my $ofile=$alnfile.'_mappedvna.fst';
$aln->export($ofile, 'fst', 1);
printf "\n- mapped folding info writing to file %s\n", $ofile;