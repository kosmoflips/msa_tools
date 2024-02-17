#!/usr/bin/perl
use strict;
use warnings;

use Data::Dumper;
use Storable qw/:DEFAULT nstore dclone/;
use File::Spec;
use Getopt::Long;

# input bed/graph files, where 1st column has gene id AND coordinates are genomic (no conversion)
# output chr and move gene id to 'name field'


my (@bed_files,$file_db,$bedgraph,$help);
GetOptions(
	"beds=s{1,}"=>\@bed_files,
	"db=s{1}"=>\$file_db,
	# "graph"=>\$bedgraph,
	"help"=>\$help,
);

if ($help or !$file_db or !@bed_files) {die <<USAGE;
-----------------------------------------
convert gene ID to chromosome in BED file

----- required -----
[-b BED_file1 2 ...] # paths to BED files to be converted
  # input bed file must have 1st column to be gene ID, and coordinates should be genomic as there is no coordinate conversion
[-d Exon_path] # this file can be:
  - a Perl-Hash file "*.hash" made by "mk_exon_from_gtf.pl"
  # SQL mode isn't currently supported!

NOTE: make sure bed file and DB file are based on the same GTF reference.
-----------------------------------------

USAGE
}


my $exondata=retrieve($file_db);
printf "\n\nloaded exon info file %s . . .\n\n", $file_db;


foreach my $file_bed (@bed_files) {
	printf "> %s . . .\n", $file_bed;
	if (!-e $file_bed or !-T $file_bed) {
		print STDERR "  file doesn't exist or isn't plain-text, skip!\n";
		next;
	}
	open (my $fh, $file_bed);

	my $ofile=$file_bed.'.remap_chr.bed';
	open (my $fh2, ">", $ofile);

	while (<$fh>) {
		if (/^track/) { # likely 1st line
			print $fh2 $_;
			next;
		}
		chomp;
		my @c=split /\t/;
		my $chr=$exondata->{$c[0]}{info}{chr}||0;
		$c[0]=($chr eq '0'? '#':'').$chr;
		printf $fh2 "%s\n", join ("\t", @c);
	}
	print "  done\n";
}