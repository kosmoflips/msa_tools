#!/usr/bin/perl
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


# parse gtf file and convert info to sqlite DB. only keep certain info (see "create tables")

my ($gtffile,$force,$resume,$help);
GetOptions(
	"gtf=s{1}"=>\$gtffile,
	"resume"=>\$resume,
	"force"=>\$force,
	"help"=>\$help,
);

# if (!$help and (!$gtffile or !-e $gtffile or -z $gtffile or $gtffile!~/\.gtf$/i)) {
	# print STDERR "gtf file not given or doesn't exists! use [-g] to specify a *.gtf file.\n";
	# exit;
# }

if ($help or !$gtffile) {die <<USAGE;
-----------------------------------------
parse GTF file and extract certain info to sqlite file

-- required --
[-g path/to/genome.gtf OR /genome.gtf.sqlite]
  - must have correct extension *.gtf or *.sqlite

-- optional --
[-f] force to rebuild gtf.sqlite file
  - will overwriting existing file under the same name
[-r] resume to finish current gtf.sqlite file
  - TRY NOT TO USE IT. NOT FULLY TESTED, BUGGY!

NOTE: [-f] overwrites [-r]

-- output --
INPUT.gtf.sqlite
  - saved in the same location
-----------------------------------------

USAGE
}

my $parse2db;
my $dbfile=0;
if ($gtffile=~/\.gtf$/i) {
	print "- parsing GTF file...\n";
	$dbfile=$gtffile.'.sqlite';
	$parse2db=1;
} else {
	printf STDERR "GTF file %s must have extension .gtf or .sqlite!\n", $gtffile;
	exit;
}

if (-e $dbfile and !-z $dbfile and !$force) {
	print "  - db file already exists\n";
	if ($resume) { # flag in case last run wasn't finished.
		$parse2db=1;
	} else {
		$parse2db=0;
	}
} else {
	$parse2db=1;
}

if ($parse2db) { # mk sqlite from gtf
	GTFsupport::gtf2db($gtffile,$dbfile,$resume);
}

printf "  - Done. GTF saved as SQLite at %s\n\n", $dbfile;
