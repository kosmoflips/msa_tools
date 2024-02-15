use strict;
use warnings;

use Data::Dumper;
use Storable qw/:DEFAULT nstore dclone/;
use File::Spec;
use File::Path;
# use File::Copy;
# use File::Temp;
use Getopt::Long;


# split a gtf file by 3rd column, which is for genomic parts, e.g. exon, cds, etc.


my @infiles;

my $help;
GetOptions(
	"inputfile|gtf=s{1,}"=>\@infiles,
	"help"=>\$help
);

if ($help or !@infiles) {
die <<HELP;
--------------------------------------------------
split input GTF files by genomic parts (3rd column)

options, * are required:
*[-g INPUT_GTF1 2 3 ...] # full path preferred
  # must have a *.gtf extension

split GTF part files will be written to subdir `gtf_parts`
--------------------------------------------------

HELP
}

foreach my $infile (@infiles) {
	printf "> %s . . .\n", $infile;
	if ($infile!~/\.gtf$/i) {
		print "  input file must have a *.gtf extension, skip!\n";
	}
	if (!-e $infile or -z $infile) {
		print "  file doesn't exist or is empty, skip!\n";
		next;
	}

	my @cpath=File::Spec->splitpath($infile);
	my $infname=pop @cpath;
	$infname=~s/\.gtf$//i;
	my $odir=File::Spec->catdir(@cpath, "gtf_parts");
	mkpath $odir if !-d $odir;

	open (my $fh, $infile);
	while (<$fh>) {
		next if !/\S/;
		next if /^#/;
		my @c=split "\t";
		my $part=$c[2];
		my $ofile=File::Spec->catfile($odir, sprintf("%s__%s.gtf", $infname, $part));
		my $fh2;
		if (!-e $ofile) {
			open ($fh2, ">", $ofile);
			printf $fh2 "# sub-GTF data for %s, extracted from original file %s.gtf\n", $part, $infname;
		} else {
			open ($fh2, ">>", $ofile);
		}
		print $fh2 $_;
	}
}


