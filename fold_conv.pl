use strict;
use warnings;

use Data::Dumper;
use File::Spec;
use File::Path;
use Getopt::Long;

use lib ".";
use GenKode::Fold;
use GenKode::Commons;

#read/convert between different folding formats: vna(fst),ct,mfold constraint

# input vna files, output to individual ct files.
# vna files may be in vna format (3 line per seq), or regular fasta,etc. format(4line per seq)

my ($help, @files, $format_out, $voff,$write_error,$reorderfile);
GetOptions(
	'help'=>\$help,
	'error'=>\$write_error,
	'voff'=>\$voff,
	'file=s{1,}'=>\@files,
	'reorder=s{1}'=>\$reorderfile,
	'x=s{1}'=>\$format_out,
);


my $filelist=GenKode::Commons::get_all_input_files(\@files);

my $outformat={
	'fst' =>1,
	'vna' =>1,
	'ct' =>1,
	'mfold' =>1,
};

$format_out=lc $format_out if $format_out;

if ($help or !$filelist or (!$reorderfile and !$format_out and !$write_error) or ($format_out and !$outformat->{$format_out})) {
	die <<HELP;

-------------------------------------
#convert folding-related formats: vna, ct, fst

>> input <<
[-f FILE1 2 DIR3 4...] # absolute path, must all have the same format
# acceptable input formats: fst, ct, vienna
  - if input file format is fst, can guess-ignore sequences that do not have a following folding line, may cause errors
# if [-d] is specified, all plain-text files will be parsed. do not read subdirs

>> output <<
[-reorder NAME_LIST_FILE] #give a text file listing all names in the new order
  # ONLY works when input file is in fasta format
  # will output to fst format
[-x OUTPUT_FORMAT] # convert input file to one of "fst, ct, vna (default), mfold"
  [-voff] #optional. only works when writing a fasta file, skip all vna lines
[-e] # check for folding errors
  # use this when you don't want to convert to ct files
  # NOT compatible with writing to mfold
  # when converting to ct files, errors will always be checked
-------------------------------------

HELP
}


foreach my $file (@$filelist) {
	printf "\n>%s. . .\n", $file;
	if (!-e $file or -z $file) {
		next;
	}
	my $m=GenKode::Fold->new;
	if (!$m->import($file)) {
		print "  - unknown format, skip. .\n";
		next;
	}
	if ($reorderfile and -e $reorderfile) { #for reorder, only deal with fst files
		my $order=getorder($reorderfile);
		# print Dumper $m->{1};
		$m->reorder($order);
		# die Dumper [keys %{$m}];
		my $ofile=$m->export($file.'_reorder', 'fst', 0, $voff);
		printf "  - sequences reordered and written to: %s\n", $ofile;
		next;
	}
	if ($format_out and $outformat->{lc $format_out}) {
		my $fout=$m->export($file, $format_out, 0, $voff);
		printf "  - convert folding to: %s\n", $fout;
	}
	if ($write_error or ($format_out and ($format_out ne 'mfold' or $format_out eq 'ct') )) {
		my $eout=$m->export_error($file);
		if ($eout) {
			printf "  - errors wrote to: %s\n", $eout;
		} else {
			print "  - no errors found\n";
		}
	}
	print "\n";
}
