#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

use GenKode::HMMER qw/extract_hmmsearch/;

# extract hmm search output into tab separated format
my @files;
my $help;
GetOptions(
	'files=s{1,}'=>\@files,
	'help'=>\$help,
);

if ($help or !@files) {
	print <<MSG;
-----------------------------------
[-f] input files, must be hmmer search output
-----------------------------------
MSG
}

foreach my $f (@files) {
	printf "%s. .\n", $f;
	if (!-e $f) {
		printf "  file not found\n";
		next;
	}
	my ($lis1,$lis2)=extract_hmmsearch($f);
	my $fout=$f.'_extracted.txt';
	open (my $fh2, ">", $fout);
	unshift @$lis2, ['']; #to add one empty line for printing
	foreach my $line (@$lis1,@$lis2) {
		printf $fh2 "%s\n", join ("\t", @$line);
	}
}