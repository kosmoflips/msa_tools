use strict;
use warnings;
use File::Spec;
use File::Path;
use Getopt::Long;
use Data::Dumper;

use lib '.';
use GenKode::MSA;
use GenKode::CentralDogma;
use GenKode::Commons;

#msa manipulation
# by kiyo @ http://www.pocchong.de
# created: 2013-autumn
# change log (not from very beginning)
# 150713 reorder
# 150730 input arg, scalar to H ref
# 150910 swap b/w d/rna
# 170120 vna2ct function removed to separate script. too complicated and requires too many flags
# 190524 simplify with huge revision (?)

my $help;
my (@files);
my (@fromformats, $toformat); #using "from", only get specified extension. omit to read all known files.
my (@mask,$sortfile,$renamefile);
my ($strip,$strip_allnogap, $usefullname);
my ($as_dna, $as_rna);

GetOptions(
	'h'=>\$help,
	'rna'=>\$as_rna,
	'dna'=>\$as_dna,
	'full'=>\$usefullname,
	'files=s{1,}'=>\@files,
	't=s{1,}'=>\@fromformats, #require word match. *.fst can't match *fsa can't match *.fasta
	'x=s{1}'=>\$toformat,
	#sort > rename > mask. no specific reason
	'mask=s{1,}'=>\@mask, #mask name, char 1, 2, ...
	'sort|reorder=s{1}'=>\$sortfile,
	'n|rename=s{1}'=>\$renamefile,
	'g'=>\$strip,
	'all'=>\$strip_allnogap,
);


# 1. grab files
my $filelist=GenKode::Commons::get_all_input_files(\@files,\@fromformats);


if ($help or !$filelist) {
	die <<USAGE;
------------------------------------------
*process MSA files.
*does NOT deal with VNA. use fold_conv.pl instead.

>>>> input <<<<
#required
[-file FILE1 2 DIR3 4...] #file/dir paths, NO subdir

#optional
[-t ReadFormat_x y ...] #only fetch files with given extension(s), e.g. -r fst phy

>>>> manipulation (optional, choose only one) <<<<
[-dna] # write all sequence as dna: U to T
[-rna] # write all sequence as dna: T to U
# will do for all sequences, so have a back up of your annotation line
# output will be fasta

[-s sort_list_path]
  - a txt file containing a list of seq-names, in desired order
  - NOTE: in input file contains a folding vna line, use fold_conv.pl instead
[-n rename_path]
  - a txt file in format of : <CURRENT-name><TAB><NEW-name>
[-g (-all)]
  - strip gaps
  - if [-all] is specified, will only retain columns that all seqs have a non-gap char
  - if [-all] is NOT specified, will retain columns that at least one seq has a non-gap char
  - IMPORTANT: ambiguous chars "N" and "X" for nt and aa are NOT treated as gaps here!
[-m MASK_Name (Mask_Char1 2 ...)] #keep columns according to sequence "MASK_Name".
  - optional: mask_char. when given will only extract specified letter/numbers in "mask" seq, can input mutiple
  - force write to fst

# NOTES:
  - [-dna] > [-rna] > [-s] > [-n] > [-g] > [-m] if more than one is given
  - [-s], [-n], [-m] are all case insensitive
  - [-s] and [-n] write to [-x] specified format or fst
  - [-g] and [-m] only write to fst

>>>> output format <<<<
[-x output_format (-full)] #use ONLY ONE of below
  - standard: fst , aln , phy, sto, nxs, mase
  - special: namelist (txt), pure (fst), paml (txt), single (fst, lose all gaps)
  - if "-full" is given, will use full name in PHY, STO block instead of 10 letters per name

------------------------------------------
USAGE
}


# 2. loop each file, do following in order, if specified
my $order;
foreach my $file (@$filelist) {
	my $m=GenKode::MSA->new;
	printf "\n>>%s . .\n", $file;
	my $readchk=$m->import($file);
	if (!$readchk) {
		print "  can't parse file. skip.\n";
		next;
	}

#---------------------------------------------

	# if given, change NT and exit
	if ($as_dna or $as_rna) {
		my $asnna=GenKode::CentralDogma->new;
		if ($as_dna) {
			$m->to_dna();
		}
		else {
			$m->to_rna();
		}
		my $ofile=$m->export('fst', {
				fileroot=>$file.($as_dna?'_asDNA.fst':'_asRNA.fst'),
				fixfilename=>1,
			});
		printf "   converted as %sNA sequence in file: %s\n", ($as_dna?'D':'R'), $ofile;
	}

#---------------------------------------------

	$toformat=$m->chk_format($toformat); #FIRST step. make sure writing is okay
	my $feval=$m->format_type($toformat);
	my $changed;
# REGARDLESS if given file works,  reorder > sort > mask
	if (!$changed and $sortfile) {
		$order=getorder($sortfile) if !$order;
		my $dosort=$m->reorder($order);
		if ($dosort) {
			$toformat='fst' if !$toformat;
			my $ofile=sprintf '%s.reorder.%s', $file, ($feval==1?$toformat:'fst');
			printf "  - reorder as in: %s\n", $sortfile;
			$changed=$ofile;
		}
	}
	if (!$changed and $renamefile) {
		my $renamemap=getrename($renamefile);
		my $rename=$m->rename($renamemap);
		if ($rename) {
			$toformat='fst' if !$toformat;
			my $ofile=sprintf '%s.rename.%s', $file, ($feval==1?$toformat:'fst');
			printf "  - rename as in: %s\n", $renamefile;
			$changed=$ofile;
		}
	}

	if (!$changed) {
		if ($strip) {
			my $outf=$m->export_stripgaps($file,$strip_allnogap);
			printf "  - output for nogap area: %s\n", $outf; #no $changed assign, as output is done here
		}
		elsif (!$changed and @mask) {
			my $outf=$m->export_mask($file,\@mask);
			printf "  - output for mask: %s\n", $outf; #no $changed assign, as output is done here
		}
	}

#---------------------------------------------

# 3. export
	#output special format if needed
	if ($feval>1) {
		if ($feval<=3) {
			my $ofile=$m->export_special($file, $toformat);
			printf "  - output for %s: %s\n", $toformat,$ofile;
		}
		elsif ($feval==4) {
			my $odir=$m->export_single($file);
			printf "  - output single seqs in: %s\n", $odir;
		}
	}
	# msa changed, output to std format anyway // otherwise only convert format if output fmt isn't the same as input file
	if ($changed or
		($feval==1 and 
				(($toformat ne $m->format)) or ($toformat eq 'fst')
			)
	) {
		my $ofile=$m->export($toformat, {
				fileroot=>$changed||$file,
				fixfilename=>$changed?1:0,
				fullname=>$usefullname,
			});
		# if ($changed) {
			# $ofile=$m->export($changed,$toformat,1);
		# } else {
			# $ofile=$m->export($file,$toformat);
		# }
		printf "  - convert to %s: %s\n", $toformat, $ofile;
		if ($toformat=~/nxs/i) {
			print "  - !!output to nxs, you may want to check molecular type, currently using 'standard'\n";
		}
	} else {
		print "  >>> done\n";
	}
}

# ------ subs --------------
sub getrename {
	my ($renamefile)=@_; #curr-name <tab> new-name
	open (my $fh, $renamefile) or return undef;
	my $map;
	while (<$fh>) {
		next if !/\S/;
		next if /^#/;
		chomp;
		my (@c)=split /\t/;
		$map->{$c[0]}=$c[1];
	}
	return $map;
}

sub verify_outputformat {
	my $toformat=shift;
	if (!$toformat) {
		$toformat='fst';
	} else {
		$toformat=lc $toformat;
	}
}

