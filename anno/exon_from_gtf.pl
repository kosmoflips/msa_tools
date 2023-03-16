#!/usr/bin/perl
use strict;
use warnings;

use Data::Dumper;
use Storable qw/:DEFAULT nstore dclone/;
use File::Spec;
use Getopt::Long;

# parse gtf file, extract "exon" records and store position info to a Perl hash

my (@gtffiles,$force,$help);
GetOptions(
	"gtf=s{1,}"=>\@gtffiles,
	"force"=>\$force,
	"help"=>\$help,
);

if ($help or !@gtffiles) {
	die <<USAGE;
-----------------------------------------
parse GTF files and save exon's position info into perl hash file

--- INPUT ---
*[-g GTF1 GTF2 ...] # one or more GTF file paths
[-f] # if output file exists, force overwriting

--- OUTPUT ---
GTF.exonpos.hash

-----------------------------------------

USAGE
}


foreach my $infile (@gtffiles) {
	printf "> %s . . .\n", $infile;
	my $ofile= $infile.'.exonpos.hash';

	if (-e $ofile and !$force) {
		printf "  - output file exist, skipped! (or use [-f] to force overwriting.)\n";
		next;
	}

	my $data;
	open (my $fh, $infile);
	my $curr_tid='';
	my $curr_trx=[];
	my $curr_chr='';
	my $curr_attr;
	while (<$fh>) {
		next if /^#/;
		chomp;
		my @c=split /\t/;
		if ($c[0] ne $curr_chr) {
			printf "  - chromosome %s, line %s . . .\n", $c[0], $.;
			$curr_chr=$c[0];
		}
		next if $c[2] ne 'exon';
		my $attr=parse_gtf_attr($c[8]);

		# gene version and chr info
		if (!exists $data->{$attr->{gene_id}}{info}) {
			$data->{$attr->{gene_id}}{info} = {
				'chr' => $c[0],
				gene_version=>$attr->{gene_version},
				strand=>$c[6] eq '-'? 2:1,
			};
		}
		
		# see if can reset trx hash
		if ($curr_tid and ($curr_tid ne $attr->{transcript_id}) ) { # has read a prev line AND current id isn't the same as prev one => new trx id
			$curr_trx=conv_exon_position($curr_trx);
			# save to master data
			$data->{$attr->{gene_id}}{$attr->{transcript_id}}=dclone $curr_trx;
			# reset
			$curr_trx=undef;
			$curr_trx=[];
			$curr_attr=dclone $attr; # copy current attr info (for the last record)
		}
		$curr_tid=$attr->{transcript_id};
		# add trx versions as [0]
		$curr_trx->[0]=$attr->{transcript_version};
		# add current exon info
		$curr_trx->[$attr->{exon_number}]=[ $c[3],$c[4] ];
	}
	# process the last id
	if ($curr_trx and $curr_attr) {
		$curr_trx=conv_exon_position($curr_trx);
		# save to master data
		$data->{$curr_attr->{gene_id}}{$curr_attr->{transcript_id}}=dclone $curr_trx;
	}
	nstore($data, $ofile);
	printf "\n  >> parsed exon position info saved to %s\n", $ofile;
}


sub parse_gtf_attr {
	# need gene/trx id and exon number
	my ($attrline)=@_;
	#gene_id "ENSG00000284662"; gene_version "1"; transcript_id "ENST00000332831"; transcript_version "4"; exon_number "1"; gene_name "OR4F16"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "OR4F16-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS41221"; exon_id "ENSE00002324228"; exon_version "3"; tag "basic"; transcript_support_level "NA (assigned to previous version 3)"
	my $attr;
	for my $item (qw/gene_id  gene_version  transcript_id  transcript_version  exon_number/) {
		if ($attrline=~/$item "(.+?)"/) {
			$attr->{$item}=$1;
		}
	}
	return $attr;
}

sub conv_exon_position {
	my ($trx)=@_;
	my $m=1;
	my $n;
	# no need to worry about strand. as long as they're properly numbered in GTF
	for my $i (1..(scalar(@$trx)-1)) { # [0] is trx id
		my ($p,$q)=@{$trx->[$i]};
		my $diff=$q-$p;
		$n=$m+$diff;
		# update current data
		push @{$trx->[$i]}, $m, $n;
		# printf "%s %s -- %s %s", $p, $q, $m, $n;<>;
		$m=$n+1;
	}
	return $trx;
}
