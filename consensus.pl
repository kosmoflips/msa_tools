use strict;
use warnings;
use File::Spec;
use File::Path;
use Getopt::Long;
use Data::Dumper;

use lib '.';
use GenKode::Fold;
use GenKode::MSA;
use GenKode::Commons;


my $help;
my (@files);
my @cutoffs;
my @exclude_file;
my $mode;

GetOptions(
	'h'=>\$help,
	'files=s{1,}'=>\@files,
	'c=f{1,}'=>\@cutoffs,
	'exclude=s{1,}'=>\@exclude_file,
	'mode=i{1}'=>\$mode, #0,1,2
);

my $filelist=GenKode::Commons::get_all_input_files(\@files);

if ($help or !$filelist) {
	die <<USAGE;
------------------------------------------
*generate consensus sequence

>>>> input <<<<
#required
[-f FILE1 2 DIR1 2...] #file or dir paths, do NOT proceed sub-dir

>>>> optional <<<<
[-m MODE] #0: only ATCGN (default), 1: with R Y, 2: all symbols
[-e EXCLUDE] # either names or ONE text file containing names to be excluded from calculating consensus, one name per line. e.g. "annotation"
[-c cutoff_value1 2 3...] between 0~1
# will always make from 0.7 to 1, 0.05 step

>>>> output <<<<
# fasta format only

------------------------------------------
USAGE
}

my $exclude;
if (@exclude_file) {
	if (-e $exclude_file[0]) {
		open (my $fh, $exclude_file[0]);
		while (<$fh>) {
			next if /^#/;
			next if !/\S/;
			chomp;
			$exclude->{lc $_}=1;
		}
	} else {
		$exclude->{lc $_}=1 foreach (@exclude_file);
		# die  Dumper $exclude;
	}
}

$mode=0 if (!$mode or $mode>2);

my $cutofflist;
foreach my $c (@cutoffs) {
	next if $c!~/^\d\.?\d*$/;
	$cutofflist->{$c}=1;
}
for (my $v=70; $v<=100; $v+=5) {
	$cutofflist->{($v/100)}=1;
}


foreach my $file (@$filelist) {
	my $m=GenKode::Fold->new;
	printf "\n- %s. .\n", $file;
	my $readchk=$m->import($file);
	# die Dumper [keys %$m];
	if (!$readchk) {
		print "  can't parse file. skip.\n";
		next;
	}
	# GenKode::Fold::parse_fold($m); #remove all vna info

	foreach my $val (sort {$a<=>$b} keys %$cutofflist) {
		$m->mk_consensus($val,$mode,$exclude);
		$exclude->{lc $m->getname($m->lastid)}=1;
	}
	my $fout=$file.'_consensus.fst';
	$m->export($fout,'fst',1);
	printf "  %s. .\n", $fout;
}
