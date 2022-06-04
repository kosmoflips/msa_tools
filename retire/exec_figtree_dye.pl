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
use GenKode::Method_Phylogeny_FigTree;

# import FigTree saved nxs file, specify groups of taxon (and colour), dye each taxa node
# WILL STILL HAVE TO check the parental node and dye manually, as this isn't easy to handle


my ($usr,$help);
GetOptions(
	"file=s"=>\$usr->{file},
	"chart=s"=>\$usr->{chart},
	"help"=>\$help,
);

my ($group,$chart)=GenKode::Method_Phylogeny_FigTree->arrange_dye($usr->{chart});

if ($help
	or (!$usr->{file} or !-e $usr->{file} or -z $usr->{file})
	or (!$usr->{chart} or !-e $usr->{chart} or -z $usr->{chart})
	or !($group and $chart)
) {
die <<USAGE;
-----------------------------------------
[-f] input file saved by FigTree
[-c] chart file for group colouring

*chart file format:
>ffcccc
taxa1 taxa2 taxa3...
>ff0066
taxa4 taxa5...
...

-----------------------------------------

USAGE
}

$usr->{fileout}=$usr->{file}.'_out.txt';
# $usr->{log}=$usr->{file}.'_log.txt';

open (my $fh, $usr->{file});
open (my $fh2,">", $usr->{fileout});
my $tf;
while (<$fh>) {
	chomp;
	if ($tf and /end;/) {
		$tf=0;
	}
	elsif (/begin trees;/) {
		$tf=1;
	}
	elsif ($tf) {
		my ($tag,$tree)=$_=~/^(\s*tree.+?=)(.+)$/;
		my $tree2=GenKode::Method_Phylogeny_FigTree->dye_tree($tree,$group,$chart);
		printf $fh2 "%s%s\n",$tag, $tree2;
		next;
	}
	printf $fh2 "%s\n", $_;
}

=pod
{
open (my $fh3,">", $usr->{log});
my $time=localtime(time);
printf $fh3 "### figtree auto-process ###
#
# %s
#
# input file: %s
# output file: %s
#
#######
",
$time, $usr->{file}, $usr->{fileout};

GenKode::Method_Phylogeny_FigTree->prt_dye_chart($group,$chart,$fh3);
}
=cut

# my $chart

