package GenKode::CentralDogma;
use strict;
use warnings;

### by kiyo @ http://www.pocchong.de

#for convenience, uc all AA, lc all NT

# 200518 : everything requires calling a "->new" method
# last update: 181207


sub new {
	my $class=shift;
	return bless {}, $class;
}

sub translate {#trans NT to AA by the given frame. if no frame given, use +1
#input nt seq should have no invalid chars/spaces/gaps etc
	# my (@things)
	my ($self, $nt,$rf,$codontable)=@_; #scarlar, digit, digit. rf is one from: +1,+2,+3,-1,-2,-3
	$rf=1 if (!$rf or $rf!~/^[+-]?[123]$/);
	$nt=_cleanseq($nt);
	$nt=$self->revcom($nt,1) if $rf<0;
	$nt=$self->d2rna($nt) if $nt=~/t/i;
	my $i=(abs $rf) + 2; #so each time the cut-off codon will have 3.
	my $aa='';
	my $codonchart=$self->codon_table($codontable);
	while ($i<=length $nt) {
		my $codon=uc (substr $nt, $i-3, 3);
		$aa.= $codonchart->{$codon}||'X';
		$i+=3;
	}
	return $aa;
}

sub revcom { #get reverse complimentary fo a given seq. return DNA unless toggled
	my ($self,$nt,$rna)=@_;
	$nt=r2dna($nt) if $nt=~/u/i;
	$nt=reverse $nt;
	$nt=~tr/atcgATCG/tagcTAGC/;
	$nt=$self->d2rna($nt) if $rna;
	return $nt;
}
sub r2dna {
	shift;
	my $nt=shift;
	$nt=~tr/uU/tT/;
	return $nt;
}
sub d2rna {
	shift;
	my $nt=shift;
	$nt=~tr/tT/uU/;
	return $nt;
}

sub find_all_orf { #input sequence, find all orf, return hash
	my ($self,$genome)=@_;
	my $start={
		'atg'=>1,
		'ttg'=>1,
		'ctg'=>1,
		'gtg'=>1,
	};
	my $glen=length $genome;
	my $orfs;
	#strand=1
	foreach my $rf (1..3) {
		my $pos=$rf-1; #initial
		my $first_min=0;
		my $first_all;
		while (($pos+3)<=$glen) {
			my $cod=substr $genome, $pos, 3;
			$cod=~s/u/t/ig;
			$cod=lc $cod;
			if ($start->{$cod}) {
				if (!$first_min) {
					$first_min=$pos+1;
				}
				if (!$first_all->{$cod}) {
					$first_all->{$cod}=$pos+1;
				}
			}
			elsif ($cod=~/^(tga|taa|tag)$/i and $first_min) {
				$orfs->{$rf}{$first_min}={
					start=>$first_min,
					end=>($pos+3),
				};
				foreach my $c (keys %$first_all) {
					$orfs->{$rf}{$first_min}{$c}=$first_all->{$c};
				}
				$first_min=0;
				$first_all=undef;
			}
			$pos+=3;
		}
	}
	#strand=2
	foreach my $rf2 (1..3) {
		my $endpos=$glen-$rf2+1-3; #initial
		my $first_min=0;
		my $first_all;
		while (($endpos-3)>=0) {
			my $cod=substr $genome, $endpos, 3;
			$cod=$self->revcom($cod);
			$cod=~s/u/t/ig;
			$cod=lc $cod;
			if ($start->{$cod}) {
				if (!$first_min) {
					$first_min=$endpos+3;
				}
				if (!$first_all->{$cod}) {
					$first_all->{$cod}=$endpos+3;
				}
			}
			elsif ($cod=~/^(tga|taa|tag)$/i and $first_min) {
				$orfs->{$rf2*-1}{$first_min}={
					start=>$first_min,
					end=>($endpos+1),
				};
				foreach my $c (keys %$first_all) {
					$orfs->{$rf2*-1}{$first_min}{$c}=$first_all->{$c};
				}
				$first_min=0;
				$first_all=undef;
			}
			$endpos-=3;
		}
	}
	return $orfs;
}

sub align_aa2nt { #input both AA and NT. if they match, align AA to NT so AA can have aligned gaps. ONE PAIR only
	my ($self, $nt, $aa, $rf, $codon_table)=@_;
	my $nt2=$nt;
	if (!$rf or $rf!~/^[+-]?[1-3]$/) {
		$rf=1;
	}
	$nt2=_cleanseq($nt);
	$aa=_cleanseq($aa,1); #in case input aa has gaps
	my $aa2=$self->translate($nt2, $rf, $codon_table);
	if (lc $aa2 eq lc $aa) { # must match. processing partial match is more complicated so I'll ignore here
		my $currpos=0;
		my $codoncount=0;
		my $aapos=0;
		my $aa_aln='';
		if ($rf>0) {
			for (my $i=0; $i<length $nt; $i++) {
				my $c=substr $nt, $i, 1;
				if ($c=~/[a-z]/i) {
					$currpos++;
					if ($currpos>=$rf) {
						$codoncount++;
						if ($codoncount==1) {
							$aa_aln.=substr $aa, $aapos, 1;
							$aapos++;
						} else {
							if ($codoncount==3) {
								$codoncount=0;
							}
							$aa_aln.='-';
						}
					} else {
						$aa_aln.='-';
					}
				} else {
					$aa_aln.='-';
				}
			}
		}
		else {
			for (my $i=(length $nt)-1; $i>=0; $i--) {
				my $c=substr $nt, $i, 1;
				if ($c=~/[a-z]/i) {
					$currpos++;
					if ($currpos>=abs($rf)) {
						$codoncount++;
						if ($codoncount==1) {
							$aa_aln.=substr $aa, $aapos, 1;
							$aapos++;
						} else {
							if ($codoncount==3) {
								$codoncount=0;
							}
							$aa_aln.='-';
						}
					} else {
						$aa_aln.='-';
					}
				} else {
					$aa_aln.='-';
				}
			}
			$aa_aln=reverse($aa_aln);
		}
		return $aa_aln;
	} else {
		return $aa;
	}
}

sub _cleanseq {
	my $seq=shift;
	my $is_aa=shift;
	if ($is_aa) {
		$seq=~s/[^a-z\*]//gi;
	} else {
		$seq=~s/[^a-z]//gi;
	}
	return $seq;
}

sub amino_abbr { #1-letter to 3 letter amino acid
	shift;
return {
A=>'Ala',R=>'Arg',N=>'Asn',D=>'Asp',
C=>'Cys',Q=>'Gln',E=>'Glu',G=>'Gly',
H=>'His',I=>'Ile',L=>'Leu',K=>'Lys',
M=>'Met',F=>'Phe',P=>'Pro',S=>'Ser',
T=>'Thr',W=>'Trp',Y=>'Tyr',V=>'Val',
};
}

# codons 
sub codon_table {#use Numbers as switch to diff codon sets. null = standard. RNA letters, UPPER CASE
# switches same as NCBI : http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes
	shift;
	my $s=shift;
	my $tb={ #std table
GCU=>'A',GCC=>'A',GCA=>'A',GCG=>'A',
CGU=>'R',CGC=>'R',CGA=>'R',CGG=>'R',AGA=>'R',AGG=>'R',
AAU=>'N',AAC=>'N',
GAU=>'D',GAC=>'D',
UGU=>'C',UGC=>'C',
CAA=>'Q',CAG=>'Q',
GAA=>'E',GAG=>'E',
GGU=>'G',GGC=>'G',GGA=>'G',GGG=>'G',
CAU=>'H',CAC=>'H',
AUU=>'I',AUC=>'I',AUA=>'I',
UUA=>'L',UUG=>'L',CUU=>'L',CUC=>'L',CUA=>'L',CUG=>'L',
AAA=>'K',AAG=>'K',
AUG=>'M',
UUU=>'F',UUC=>'F',
CCU=>'P',CCC=>'P',CCA=>'P',CCG=>'P',
UCU=>'S',UCC=>'S',UCA=>'S',UCG=>'S',AGU=>'S',AGC=>'S',
ACU=>'T',ACC=>'T',ACA=>'T',ACG=>'T',
UGG=>'W',
UAU=>'Y',UAC=>'Y',
GUU=>'V',GUC=>'V',GUA=>'V',GUG=>'V',
UAA=>'*',UGA=>'*',UAG=>'*',
};
	if (!$s or $s!~/^\d+$/) { ;}
	elsif ($s==2) { # The Vertebrate Mitochondrial Code
		$tb->{AGA}='*';
		$tb->{AGG}='*';
		$tb->{AUA}='M';
		$tb->{UGA}='W';
	}
	elsif ($s==3) { # Yeast Mitochondrial Code
		$tb->{AUA}='M';
		$tb->{CUU}='T';
		$tb->{CUC}='T';
		$tb->{CUA}='T';
		$tb->{CUG}='T';
		$tb->{UGA}='W';
		$tb->{CGA}=undef;
		$tb->{CGC}=undef;
	}
	elsif ($s==4) {#The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
		$tb->{UGA}='W';
	}
	elsif ($s==5) {# Invertebrate Mitochondrial Code
		$tb->{AGA}='S';
		$tb->{AGG}='S';
		$tb->{AUA}='M';
		$tb->{UGA}='W';
	}
	elsif ($s==6) {#Ciliate, Dasycladacean and Hexamita Nuclear Code
		$tb->{UAA}='Q';
		$tb->{UAG}='Q';
	}
	elsif ($s==9) {#Echinoderm and Flatworm Mitochondrial Code
		$tb->{AAA}='N';
		$tb->{AGA}='S';
		$tb->{AGG}='S';
		$tb->{UGA}='W';
	}
	elsif ($s==10) {#The Euplotid Nuclear Code
		$tb->{UGA}='C';
	}
	elsif ($s==11) {#the Bacterial, Archaeal and Plant Plastid Code
		;#codon itself is same as No.1
	}
	elsif ($s==12) {#The Alternative Yeast Nuclear Code
		$tb->{CUG}='S';
	}
	elsif ($s==13) {#The Ascidian Mitochondrial Code
		$tb->{AGA}='G';
		$tb->{AGG}='G';
		$tb->{AUA}='M';
		$tb->{UGA}='W';
	}
	elsif ($s==14) {#The Alternative Flatworm Mitochondrial Code
		$tb->{AAA}='N';
		$tb->{AGA}='S';
		$tb->{AGG}='S';
		$tb->{UAA}='Y';
		$tb->{UGA}='W';
	}
	elsif ($s==15) {#Chlorophycean Mitochondrial Code
		$tb->{TAG}='L';
	}
	elsif ($s==21) {# Trematode Mitochondrial Code
		$tb->{TGA}='W';
		$tb->{ATA}='M';
		$tb->{AGA}='S';
		$tb->{AGG}='S';
		$tb->{AAA}='N';
	}
	elsif ($s==22) {#Scenedesmus obliquus Mitochondrial Code
		$tb->{TCA}='*';
		$tb->{TAG}='L';
	}
	elsif ($s==23) {#Thraustochytrium Mitochondrial Code 
		$tb->{UUA}='*';
	}
	elsif ($s==24) {#Pterobranchia Mitochondrial Code
		$tb->{AGA}='S';
		$tb->{AGG}='K';
		$tb->{UGA}='W';
	}
	elsif ($s==25) {#Candidate Division SR1 and Gracilibacteria Code
		$tb->{UGA}='G';
	}
	return $tb;
}
1;
