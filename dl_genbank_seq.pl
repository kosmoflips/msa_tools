#!/usr/bin/perl
use strict;
use warnings;

use File::Spec;
use File::Copy;
use File::Path;
use Getopt::Long;

use lib '.';
use GenKode::GenBank;
# use GenKode::Method_MSA;

### ref link: http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
## get gi from accession num: http://www.ncbi.nlm.nih.gov/nuccore/ABIY02000120.1?report=gilist&log$=seqview&format=text
## sample fetch url: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=189430232&strand=1&seq_start=1&seq_stop=600&rettype=fasta&retmode=text');


my ($csv,$odir,$help,$type,$mode,$db,$merge,$force);
GetOptions(
	'list=s'=>\$csv,
	'help'=>\$help,
	'odir=s'=>\$odir,
	'type=s'=>\$type,
	'mode=s'=>\$mode,
	'db=s'=>\$db,
	'merge'=>\$merge,
	'force'=>\$force,
);
if ($help or !$odir or !$csv or (!-e $csv or -z $csv or !-r $csv)) { die <<USAGE;
-------------------------
[-odir OUTPUT_DIR] #save downloads to DIR...
[-db] #output DB: nuccore (default), protein, ...
[-type] #output TYPE: fasta (default), gb, asn, ...
[-mode] #output MODE: xml, txt (default), ...
[-merge] #if given, will combine all downloaded seqs to one file (not recommend if getting whole genomes)
[-force] #if given, will overwrite existing files if the script assigns a same file name

[-list CSV_PATH] #download data in csv. FORMAT CORRECTLY!
==========
  use TAB as separator
  NO header
  follow this column order:
 >> <ACCESSION>(<START><END><STRAND><NAME>)
   - ACCESSION or GI. Accession should have version, otherwise will use genbank's current default
   - STRAND: 1 => plus, -1 or 2 => minus. default=1
   - NAME: usr-def to name the downloaded file

 >> COMMON DB/TYPE/MODE
  DB: nuccore, nucest, nucgss, protein, popset
  output format:
       GenBank XML: TYPE => gb, MODE => xml (*do NOT use for DB=>protein)
       GenBank flat: TYPE => gb, MODE => text (*do NOT use for DB=>protein)
       FASTA: TYPE => fasta, MODE => text (DEFAULT)

#more info:
    http://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.chapter4_table1/?report=objectonly
==========

#written only for download nt/pr seq/entries. do NOT use for others like PubMed!!
-------------------------

USAGE
}

mkpath $odir unless -d $odir;
my $errfile=File::Spec->catfile($odir,"_errors.txt");
open (my $ERRFH,">", $errfile);

my $mergefile=File::Spec->catfile($odir,"_MergedChunk.txt");
my $MERFH;
if ($merge) {
	open ($MERFH,">", $mergefile);
}

$db='nuccore';

my $i=0;my $e=0;
my $countall=0;
open (my $fh, $csv) or die "can't open csv file!";
while (<$fh>) { #b/c sometimes the input csv may have a crazy size
	chomp;
	next if /^#/;
	$i++;
	my @raw=split /\t/,$_;
	printf "%d. %s - ",$i , $raw[0];
	#fetch the file by request url
	my $saved=fetchseq({
		acc=>$raw[0], # accession
		start=>($raw[1]?$raw[1]:''),
		end=>($raw[2]?$raw[2]:''),
		strand=>($raw[3]?$raw[3]:''),
		name=>($raw[4]?$raw[4]:''),
		db=>$db,
		type=>($type?$type:'fasta'),
		mode=>($mode?$mode:'txt'),
		odir=>$odir,
	}, $force);
	if (!$saved or !-e $saved or -z $saved) {
		printf $ERRFH "%s\n", $_;
		printf "%s\n", $_;
		$e++;
		next;
	} else {
		printf "success\n";
	}
	#merge to master file
	if ($merge) {
		open (my $fh, $saved);
		while (<$fh>) {
			if (/^>/) {
				chomp;
				printf $MERFH ">%s\n", ($raw[4]||$');
			} else {
				print $MERFH $_;
			}
		} #copy the file to master file
		print $MERFH "\n"; #separate files
	}

	#move the file
	# my @fname=File::Spec->splitpath($saved);
	# move ($saved,File::Spec->catfile($odir,pop @fname));

	$countall++;
	sleep 2; #reduce server load
	if ($countall%100==0) {
		print "sleep 30 seconds to reduce remote server load...\n\n";
		sleep 30;
	}
}
printf "\nTotal sequence: %i\nFailed: %i",$i,$e;
close ($ERRFH);
unlink $errfile if -z $errfile;
