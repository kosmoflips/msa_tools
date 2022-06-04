use strict;
use warnings;

use Data::Dumper;
use File::Spec;
use File::Path;
use Getopt::Long;

use lib ".";
use GenKode::Fold;

# input MSA files and folding info files (vna,ct), add folding to MSA
# must have exactly the same sequence title, and sequence must match exactly (case follow MSA)
# output to only fasta format, do NOT check folding error

my ($help, @files, @folds,);
GetOptions(
	'help'=>\$help,
	'input|file=s{1,}'=>\@files,
	'fold=s{1,}'=>\@folds,
);

my $outformat={
	'fst' =>1,
	'vna' =>1,
	'ct' =>1,
	'mfold' =>1,
};

my $m=GenKode::Fold->new;
my $m2=GenKode::Fold->new;
$m->import('D:\kosmoflips\tmp\testfold.fst');
$m2->import('D:\kosmoflips\tmp\t010_c039_mfold01.ct');

$m->match_fold(1,GenKode::Fold::mkfold('gggaaaccc','(((...)))'));


__END__

if ($help or !@files or (!$outformat->{lc $format_out} and !$write_error)) {
	die <<HELP;

-------------------------------------
#convert folding-related formats: vna, ct, fst

>> input <<
[-f FILE1 FILE2 ...] # absolute path, must all have the same format
[-d DIR1 DIR2 ...] # absolute path, will pick up all files in specified input_format (good for a bunch of .ct files)
# use at least one of [-f]/[-d] for input
# acceptable input formats: fst, ct, vienna
  - if input file format is fst, can guess-ignore sequences that do not have a following folding line, may cause errors
# if [-d] is specified, all plain-text files will be parsed. do not read subdirs

>> output <<
[-x OUTPUT_FORMAT] # convert input file to one of "fst, ct, vna (default), mfold"
  [-voff] #optional. only works when writing a fasta file, skip all vna lines
[-e] #analyse all foldings and write errors
# use at least one of [-e] and [-x]
# output to ct will always write errors
-------------------------------------

HELP
}


foreach my $file (@files) {
	printf ">%s. . .\n", $file;
	if (-e $file and !-z $file) {
		my $m=GenKode::Fold->new;
		if (!$m->import($file)) {
			print "  - unknown format, skip. .\n";
			next;
		}

		if ($format_out and $outformat->{lc $format_out}) {
			my $fout=$m->export($file, $format_out, 0, $voff);
			printf "  - convert folding to: %s\n", $fout;
		}
		if (($write_error and $format_out!~/mfold/i) or $format_out=~/ct/i) {
			my $eout=$m->export_error($file);
			if ($eout) {
				printf "  - errors wrote to: %s\n", $eout;
			} else {
				print "  - no errors found\n";
			}
		}
	}
}
