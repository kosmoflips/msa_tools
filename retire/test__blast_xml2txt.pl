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

use lib ".";
use GenKode::Commons;
use GenKode::BlastParse;

# convert xml to txt result format


my ($help,@xmlfile,$merge);
GetOptions(
	"xmlfile=s{1,}"=>\@xmlfile,
	"merge"=>\$merge,
	"help"=>\$help,
);

my $files=get_all_input_files(\@xmlfile,'xml');

if ($help or !$files) {die <<USAGE;
-----------------------------------------
# convert blast xml2 result to txt format

[-x XML_FILE1 2 ...] # must be in BLAST's "xml2" format
[-m] #optional. merge multiple xml files into one text file (when multiple xml files given). default=no merging

-----------------------------------------

USAGE
}

my $ofile1;
if ($merge) {
	$ofile1=$files->[0].'.txt';
}
foreach my $f (@$files) {
	my $ofile=$merge?$ofile1:$f.'.txt';
	open (my $fh2, ">", $ofile);
	printf $fh2 "# ------- source file: %s -------\n\n", $f;
	my $xr=GenKode::BlastParse->new($f);
	$xr->xml2txt($fh2);
}
