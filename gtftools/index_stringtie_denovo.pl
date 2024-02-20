use strict;
use warnings;

use Data::Dumper;
use Storable qw/:DEFAULT nstore dclone/;
use File::Spec;
# use File::Path;
# use File::Copy;
# use File::Temp;
use Getopt::Long;

use lib ".";
use GTFsupport;

# index denovo assembly
# input : denovo assembly generated GTF by stringtie
# output : a perl hash file, listing all start/end coordinates for each row

my @gtfs;
my $help;
GetOptions(
	"inputfile|gtf=s{1,}"=>\@gtfs,
	"help"=>\$help
);

if ($help or !@gtfs) {
die <<HELP;
--------------------------------------------------
** this script should be run from its storing directory **

index a Stringtie-generated GTF file for each exon.

----- required -----
[-g INPUT_GTF] # full path preferred
 - extension must be *.gtf

output file is a perl hash and can be read by Perl's Storable::retrieve,
access the hash by key `_gene.id_` to view the data structure.
--------------------------------------------------

HELP
}

foreach my $file (@gtfs) {
	printf ">>%s . . .\n", $file;
	if (!-e $file) {
		print "  file doesn't exist, skip!\n";
	}
	if ($file!~/\.gtf$/i) {
		print "  file must have a *.gtf extension, skip!\n";
	}
	my $idx=stringtie_gtf_indexer($file);
	my $ofile=$file.'.index.hash';
	printf "  writing indexed coordinate info to %s . . .\n", $ofile;
	nstore($idx, $ofile);
}

print "\n\nall done.";