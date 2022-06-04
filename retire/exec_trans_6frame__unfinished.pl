use strict;
use warnings;

use File::Spec;
# use File::Copy;
use Getopt::Long;

use lib ".";
use GenKode::CentralDogma;
use GenKode::Commons;
use GenKode::MSA;

use Data::Dumper;

my (@input,$table_id,@rfs,$help,$plain);
GetOptions(
	"file=s{1,}"=>\@input,
	"rfs=i{1,}"=>\@rfs,
	"plain=s"=>\$plain,
	# "min=i"=>\$minlen, #only find orfs with specifid # of AA, excluding stop
	"table=i"=>\$table_id,
	"help"=>\$help,
);

my $filelist=GenKode::Commons::get_all_input_files(\@input)||undef;

if ($help or !$filelist) {
die <<USAGE;

-----------------------------------------
translate given sequence(s) in 6 frame with coordinates

[-file] # REQUIRED. input files or dirs of target files, fasta format only

[-rf] #from -3 to 3, if not given will translate all 6 frames
[-plain]
    #if not given, will output aligned translated AA, original NT alignment will be expanded to allowing 3-letter coden space
    #if given, will only output unaligned translated AA

[-min] # only find ORFs with at least specified number of 
AA, excluding stop codon
[-table_id] # translate table id as defined on NCBI, from 1~25 non-continuous, default=1
-----------------------------------------

USAGE
}

my $rf2;
foreach my $rf (@rfs) {
	if ($rf>=-3 and $rf<=3 and $rf!=0) {
		push @$rf2, $rf;
	}
	else {
		$rf2=[qw/1 2 3 -1 -2 -3/];
	}
}

# $minlen=1 if !$minlen;
if (!$table_id or $table_id!~/^\d+$/ or $table_id>25) {
	$table_id=1;
}

my $trans=GenKode::CentralDogma->new;
foreach my $file (@$filelist) {
	printf ">%s...\n",$file;
	my $m=GenKode::MSA->new;
	$m->import($file);
	my $file2=$file.'__frame_translate.fst';
	open (my $fh2, ">", $file2);
	foreach my $id (1..$m->lastid) {
		if (!$plain) { # not plain mode, need to print both NT and AA
			printf $fh2 ">%s\n%s\n", $m->getname($id), $m->getseq($id);
		}
		my $seq=$m->getseq($id);
		foreach my $rf (@$rf2) {
			my $sa=$trans->translate($seq, $rf);
			if (!$plain) { # non plain mode, need to align translated AA
				my $sa2=$trans->align_aa2nt($seq, $sa, $rf);
				$sa=$sa2;
			}
			printf $fh2 ">%s__RF=%s\n%s\n", $m->getname($id), $rf, $sa;
		}
	}
}
