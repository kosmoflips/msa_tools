#!/usr/bin/perl
use strict;
use warnings;

use Data::Dumper;
use Storable qw/:DEFAULT nstore dclone/;
use File::Spec;
use Getopt::Long;

use lib ".";
use GTFsupport;

# parse gtf file, extract "exon" records and store position info to a Perl hash

# lighter than converting whole GTF to DB. useful if you only need to convert between transcriptomic to genomic coordinates

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

--- required ---
[-g GTF1 GTF2 ...] # one or more GTF file paths

-- optional --
[-f] # if output file exists, force overwriting

--- output ---
GTF.exons.hash

-----------------------------------------

USAGE
}


foreach my $infile (@gtffiles) {
	printf "> %s . . .\n", $infile;
	my $ofile= $infile.'.exons.hash';

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
		my $attr=GTFsupport::parse_gtf_attr($c[8]);

		# gene version and chr info
		if (!exists $data->{$attr->{gene_id}}{info}) {
			$data->{$attr->{gene_id}}{info} = {
				'chr' => $c[0],
				gene_version=>$attr->{gene_version},
				strand=>$c[6] eq '-'? 2:1,
			};
		}

		if (!$data->{$attr->{gene_id}}{$attr->{transcript_id}}[0]) { # transcript ver
			$data->{$attr->{gene_id}}{$attr->{transcript_id}}[0]={transcript_version=>$attr->{transcript_version}};
			# die Dumper $data->{$attr->{gene_id}};
		}
		$data->{$attr->{gene_id}}{$attr->{transcript_id}}[$attr->{exon_number}]=[ $c[3],$c[4] ];
	}
	foreach my $gid (keys %$data) {
		foreach my $tid (keys %{$data->{$gid}}) {
			my $tmp=$data->{$gid}{$tid};
			next if (ref $tmp)!~/array/i; # skip {info}
			$data->{$gid}{$tid}=dclone map_exon_position($tmp);
		}
	}
	nstore($data, $ofile);
	printf "\n  >> parsed exon position info saved to %s\n", $ofile;
}
