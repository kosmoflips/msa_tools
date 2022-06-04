#!/usr/bin/perl
use strict;
use warnings;

use File::Spec;
use Storable qw/:DEFAULT nstore dclone/;
use Data::Dumper;
use Getopt::Long;

use lib ".";
use GenKode::Commons;
use GenKode::BlastParse;

#extract all hit sequence (protein ID) from blast xml, first search only
#for extracting certain hit/hsp, use individual subroutine in package

my ($help,@xmlfile,$queryfile);
GetOptions(
	"xmlfile=s{1,}"=>\@xmlfile,
	"queryfile=s{1}"=>\$queryfile,
	"help"=>\$help,
);

if ($help or !@xmlfile or !chkfile($queryfile)) {die <<USAGE;
-----------------------------------------
# extract sequence alignment from given blast output xml file(s)

*required:
[-q QUERY_FILE] # either plain text or fasta
[-x BLAST_OUTPUT_FILE1 FILE2 ...] # MUST be in "xml2" format

IMPORTANT/your responsibility:
  make sure query sequence and xml files exactly match

-----------------------------------------

USAGE
}


my $query=get1stquery($queryfile);
my $ofile2=$xmlfile[0]."_printaln.fst";
open (my $fh2, ">", $ofile2) or die "can't open file to write output aln";
print "\n";
print '-' x 50, "\n";
xmlfile_extract_aln(\@xmlfile, $query, $fh2);
printf "\n\ndone. output written to file: %s\n", $ofile2;
print '-' x 50, "\n\n";

# -----subs----------------------
sub get1stquery {
	my ($f)=@_;
	if (chkfile($f)) {
		my $q='';
		open (my $fh, $f);
		while (<$fh>) {
			chomp;
			if (/^>/) {
				if ($q) {
					last;
				}
			}
			else {
				$q.=$_;
			}
		}
		$q=~s/\d|\s//g;
		return $q;
	} else {
		return '';
	}
}
sub xmlfile_extract_aln { #non-OO SHELL
	my ($xml_files, $query, $fh2)=@_;
	if (!$fh2) {
		$fh2=*STDOUT;
	}
	$query=~s/\W|_//g;
	my $qlen=length $query;

	#read all xml files
	my $xrs;
	foreach my $xf (@$xml_files) {
		my $xr2;
		printf "reading file %s. . .\n", $xf;
		eval { $xr2=GenKode::BlastParse->new($xf); };
		next if $@;
		my @fns=File::Spec->splitpath($xf);
		$xrs->{$fns[-1]}=dclone $xr2;
	}

	#process data
	my $x0=GenKode::BlastParse->new; #empty shell
	print "extracting alignment data. . .\n";
	my $alns=$x0->extract_aln($xrs);
	print "calculating column width. . .\n";
	my $colmax=$x0->calc_col_len($alns, $qlen);

	#print to FH
	$x0->print_aln_query($query,$colmax,$fh2);
	foreach my $fname (keys %$alns) {
		printf "printing %s . . .\n", $fname;
		printf $fh2 ">%s\nxxxxx\n", $fname;
		$x0->print_aln_hsp_all($alns->{$fname},$colmax,$fh2);
	}
	1;
}
