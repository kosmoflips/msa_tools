#!/usr/bin/perl
use strict;
use warnings;

use Data::Dumper;
use File::Spec;
use Getopt::Long;
use Storable qw/:DEFAULT nstore dclone/;

use lib ".";
use GTFsupport;

# input:
# - a bed file, which has its "chrom" column as transcript ID, start/end are relative to the transcript
# - a sqlite/hash file (produced by sqlite_from_gtf OR exon_from_gtf.pl)
# make sure to use exactly the same gtf since there won't be any transcript version match!

# output:
# bed file, with the "chrom" column filled with chromosome info, start/end convered to genomic, atranscript ID, other info in name field


my (@bed_files,$file_db,$bedgraph,$help);
GetOptions(
	"beds=s{1,}"=>\@bed_files,
	"db=s{1}"=>\$file_db,
	"graph"=>\$bedgraph,
	"help"=>\$help,
);

if ($help or !$file_db or !@bed_files) {die <<USAGE;
-----------------------------------------
convert transcript ID/coordinates to genomic of BED file

----- required -----
[-b BED_file1 2 ...] # paths to BED files to be converted
  # input bed file must have 1st column to be transcript ID, and coordinates should be related to transcript IDs
[-d DB/Exon_path] # this file can be:
  - a SQLite file "*.sqlite" made by "mk_sqlite_from_gtf.pl", OR
  - a Perl-Hash file "*.hash" made by "mk_exon_from_gtf.pl"

----- optional -----
[-g] # use this flag if input file is of the bedgraph format, which has only 4 columns

NOTE: make sure bed file and DB file are based on the same GTF reference.
-----------------------------------------

USAGE
}

if (-e $file_db and !-z $file_db and $file_db!~/\.sqlite|\.hash|/i) {
	print STDERR "the given <file_db> is unlikely to be a SQLite/perl-hash file. Retry.\n";
	exit;
}

my $dbmode=0;
if ($file_db!~/\.hash$/i) {
	$dbmode=1;
}

my $dbh;
my $DATA={};
my $exondata;
if ($dbmode) {
	$dbh=connectdb($file_db);
	printf "\n\nconnected to db file %s . . .\n\n", $file_db;
} else {
	$exondata=retrieve($file_db);
	printf "\n\nloaded exon info file %s . . .\n\n", $file_db;
}

foreach my $file_bed (@bed_files) {
	printf "> %s . . .\n", $file_bed;
	my $time1=time;
	printf "** start time:\%s **\n", $time1;
	if (!-e $file_bed or !-T $file_bed) {
		print STDERR "  file doesn't exist or isn't plain-text, skip!\n";
		next;
	}
	open (my $fh, $file_bed);
	my $ofile=$file_bed.'.remap_genomic.bed';
	open (my $fh2, ">", $ofile);

	while (<$fh>) {
		if ($.%100000==0) { # track time every X lines
			my $time2=time;
			printf "line %s [%.1f min] . . .\n", $., ($time2-$time1)/60;
		}
		if (/^track /) { # skip track lines, they're used to format bed files
			print $fh2 $_;
			next;
		}
		chomp;
		my @c=split /\t/;
		if (scalar @c <3) { # simple filter to remove lines that don't have at least required "chr/start/end" columns
			next;
		}

		# rm ver from id
		my ($trid,$trver)=split /\./, $c[0];
		# get start/end pos in trx
		my $p=$c[1]+1; # from 0-based to 1-based
		my $q=$c[2]; # because BED 's end position isn't included, and is 0-based. do nothing here so the $q I'm using is 1-based and indicates the true ending point
		printf "[%s] %s : %s-%s\n", $., $trid, $p, $q;
		
		# fetch pos in DB, convert to genomic
		my $exons;
		if ($dbmode) {
			$exons=get_exon_data_from_sqlite($dbh, $trid, $DATA);
		} else {
			$exons=get_exon_data_from_exons($exondata, $trid);
		}
		# print Dumper $exons,234;exit;
		my $coord=conv_coord_trx2gen($p, $q, $exons);
		if ($coord>0) { # info is converted without problem
			my $chr=$exons->[0]{chr};
			foreach my $line (@$coord) {
				# write to new bed file, start_pos converts to 0-based , end_pos is exclusive, so no change
				# not sure how BAM's score correspond to a 1000-based scale in bed, and how strand is proceeded in BAM. just keep same info as Bam2bed output
				my $pcsname=sprintf "%s:%d-%d:exon%d", $c[0], $p, $q, $line->[0];
				if ($bedgraph) { # only need 4 columns max
					printf $fh2 "%s\n", join "\t", ($chr, ($line->[1]-1), $line->[2], $c[4]||''); # 4th col is score
				} else {
					printf $fh2 "%s\n", join "\t", ($chr, ($line->[1]-1), $line->[2], $pcsname, $c[4]||'', $c[5]||''); # unsure how the strand info was from BAM but I'm just copying the score and strand cols
				}
			}
		} else { # trx id isn't found in db
			print $fh2 "##transcript_id_not_found_in_database\n";
			print $fh2 "# ", $_, "\n";
		}

		# chk tmp storage size. reset if needed
		if ($dbmode and (scalar keys %$DATA)==2500) { # only store N temp gene records in memory
			$DATA=undef;
			$DATA={};
		}
	}
	printf "** finish time: %s **\n", time;
}

