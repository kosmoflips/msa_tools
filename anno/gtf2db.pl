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

use File::Basename;
use lib "."; # so make sure gtfDB.pm is in current running dir
use mkAnnoDB;

# parse gtf file into sqlite file. only keep certain info (see "create tables")

my ($gtffile,$force,$resume,$help);
GetOptions(
	"gtf=s{1}"=>\$gtffile,
	"resume"=>\$resume,
	"force"=>\$force,
	"help"=>\$help,
);

if (!$help and (!$gtffile or !-e $gtffile or -z $gtffile or $gtffile!~/\.gtf$/i)) {
	print STDERR "gtf file not given or doesn't exists! use [-g] to specify a *.gtf file.\n";
	exit;
}

if ($help or !$gtffile) {die <<USAGE;
-----------------------------------------
parse GTF file and extract certain info to sqlite file

[-g path/to/genome.gtf OR /genome.gtf.sqlite]
*[-f] force to rebuild gtf.sqlite file. *optional*
*[-r] resume to finish current gtf.sqlite file--may be buggy! *optional*

NOTE: [-f] overwrites [-r]
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
	gtf2db($gtffile,$dbfile,$resume);
}

printf "  - Done. GTF saved as SQLite at %s\n\n", $dbfile;

