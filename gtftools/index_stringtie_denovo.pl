use strict;
use warnings;

use Data::Dumper;
use Storable qw/:DEFAULT nstore dclone/;
use File::Spec;
# use File::Path;
# use File::Copy;
# use File::Temp;
use Getopt::Long;

use lib ".";
use GTFsupport;

# index denovo assembly
# input : denovo assembly generated GTF by stringtie
# output : a perl hash file, listing all start/end coordinates for each row

my @gtfs;
my $help;
GetOptions(
	"inputfile|gtf=s{1,}"=>\@gtfs,
	"help"=>\$help
);

if ($help or !@gtfs) {
die <<HELP;
--------------------------------------------------
** this script should be run from its storing directory **

index a Stringtie-generated GTF file for each exon.

----- required -----
[-g INPUT_GTF] # full path preferred
 - extension must be *.gtf

output file is a perl hash and can be read by Perl's Storable::retrieve,
access the hash by key `_gene.id_` to view the data structure.
--------------------------------------------------

HELP
}

foreach my $file (@gtfs) {
	printf ">>%s . . .\n", $file;
	if (!-e $file) {
		print "  file doesn't exist, skip!\n";
	}
	if ($file!~/\.gtf$/i) {
		print "  file must have a *.gtf extension, skip!\n";
	}
	open (my $fh, $file);
	my $idx={
		"_gene.id_" => {
			"_ensembl_" => [ 'ref_gene_id', 'ref_gene_name', 'chr', 'strand'],
			'transcript.id' => [
				[ 'ref_trx_id' ],
				[ 'exon1_genome_start', 'genome_end', 'trx_start', 'trx_end' ],
				[ 'exon2_genome_start', 'genome_end', 'trx_start', 'trx_end' ],
			]
		}
	};
	my $cp=0;
	my $cq=0;
	while (<$fh>) {
		next if /^#/;
		next if !/\S/;
		chomp;
		my @c=split /\t/;
		next if $c[2] eq 'transcript'; # only focus on exon info in stringtie-GTF
		my $attr=parse_gtf_attr($c[-1]);
		if (!$idx->{$attr->{gene_id}}) { # this row is a new STRG.id, add gene info to key {ensembl}
			# as a A-ref: ref-gene-id and ref-gene-name
			$idx->{$attr->{gene_id}}{_ensembl_}=[ $attr->{ref_gene_id}||undef, $attr->{ref_gene_name}||undef, $c[0], ($c[6] eq "+"?1:2) ]; # ref gene id, ref gene name, strand
		}
		if ($attr->{exon_number}==1) { # this row is a new transcript, add trx info
			$idx->{$attr->{gene_id}}{$attr->{transcript_id}}[0]=[ $attr->{reference_id}||undef ]; # ref transcript id
			$cp=1; # reset trx start for new exon1
		} else {
			$cp=$cq+1; # current exon start will be +1 from previous exon's end
		}
		$cq=$cp+abs($c[3]-$c[4]);
		# add genomic coordinates for this exon, here, relate the index number in the A-ref
		$idx->{$attr->{gene_id}}{$attr->{transcript_id}}[$attr->{exon_number}]=[$c[3], $c[4], $cp, $cq]; # chr, genome-start, genome-end, trx-start, trx-end
	}

	my $ofile=$file.'.index.hash';
	printf "  writing indexed coordinate info to %s . . .\n", $ofile;
	nstore($idx, $ofile);
}

print "\n\nall done.";