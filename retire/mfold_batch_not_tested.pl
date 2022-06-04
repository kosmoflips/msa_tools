use strict;
use warnings;

use File::Spec;
use File::Path;
use Data::Dumper;
use Getopt::Long;

use lib ".";
use GenKode::Commons;
use GenKode::Fold;

#note. mfold seems hate long file path/name, so keep it as short as possible, best to chdir to reduce the length too

# input dir/file, only take fst files, for each seq, run mfold with default settings
# output will be a newly created folder in the same path

my (@files);
my $dna;
my $help;

GetOptions(
	'h'=>\$help,
	'files=s{1,}'=>\@files,
	'dna'=>\$dna,
);

my $okayinformat=[qw/fsa fst fas fasta/];
my $filelist=GenKode::Commons::get_all_input_files(\@files,$okayinformat);

if ($help or !$filelist) {
	die <<USAGE;
------------------------------------------
!!! make sure mfold is installed and can be called by typing "mfold"

[-file FILE1 2 DIR3 4...] #do NOT proceed sub-dir
# will only process fasta files. extension: *.fst *.fsa *.fas *.fasta

[-dna] fold all sequences as DNA, otherwise as RNA (default)
------------------------------------------

USAGE
}

foreach my $file (@$filelist) {
	my $m=GenKode::Fold->new;
	$m->import($file);
	printf "- %s . .\n", $file;
	my @d0=File::Spec->splitpath($file);
	my $fname=pop @d0;
	my $odir1=File::Spec->catdir(@d0, $fname.'_mfold');
	foreach my $id (1..$m->lastid) {
		printf "  - %s\n", $m->getname($id);
		open (my $fh1, ">", $tmpfile);
		print $fh1 $m->getseq($id);
		close ($fh1);
		my $odir2=File::Spec->catdir($odir1,(sprintf '%03d',$id));
		mkpath $odir2 if !-d $odir2;
		chdir $odir2;
		my $nfile=File::Spec->catfile($odir2,'info.txt');
		open (my $fh1b, ">", $nfile);
		$m->write_fh_fst_1seq($fh1b, $m->getname($id),$m->getseq($id));
		my $cmd=sprintf 'mfold SEQ=%s%s', $nfile, ($dna?' NA=DNA':'');
		eval {system ($cmd); };
	}
}

