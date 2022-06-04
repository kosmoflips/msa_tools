package GenKode::HMMER;
use strict;
use warnings;

use File::Spec;
use File::Path;
use File::Temp;
use Data::Dumper;
use Storable qw/dclone nstore/;


use base 'Exporter';
our @EXPORT = qw/

/;

our @EXPORT_OK = qw/
extract_hmmsearch
/;

our %EXPORT_TAGS = (
# translation=>[qw/ /]
);


my $HMMHEADER1=[qw/  full-ev	 full-score	 full-bias	best1-ev	best1-score	best1-bias	dom-exp	dom-Num 	seq  desc /];
my $HMMHEADER2=[qw/  seq  num  |  score  bias  c-Evalue  i-Evalue hmm-from  hmm-to  |   ali-from  ali-to  |  env-from  env-to   |  acc/];

sub extract_hmmsearch { #for hmmsearch etc. extract score table && domain annotation to table
	my $file=shift;
	if (!-e $file) {
		return undef;
	}

	open (my $fh, $file);
	my ($in, $in2);
	my $name2;
	my ($list,$list2);
	push @$list, (dclone $HMMHEADER1);
	push @$list2, (dclone $HMMHEADER2);
	unshift @$list, ["# Scores for complete sequences (score includes all domains)"];
	unshift @$list2, ["# Domain annotation for each sequence"];
	while (<$fh>) {
		chomp;
		next if /^#/;
		next if !/\S/;
		next if /^(-|\s)+$/;
		if (/Internal pipeline statistics summary/i) {
			last;
		}
		elsif (/Domain annotation for each sequence/i) {
			$in=0;
			$in2=1;
		}
		elsif (/bias/i and !$in2) {
			$in=1;
		}
		elsif ($in) {
			next if /inclusion threshold/;
			my ($line)=$_=~/^\s*(.+)\s*$/;
			my (@c)=split /\s+/, $line, 10; #split "10" times as of hmm3 format, so the last column is for description.
			push @$list, [@c];
		}
		elsif ($in2) {
			if (/Alignments for each domain/) {
				$name2=undef;
			}
			elsif (/^>>\s*/) {
				$name2=$';
				$name2=~s/\s+$//g;
			}
			elsif (/^\s*\d/ and $name2) {
				my ($line)=$_=~/^\s*(.+)\s*$/;
				my @c2=split /\s+/, $line;
				push @$list2, [$name2, @c2];
			}
		}
	}
	return ($list,$list2);
}
1;