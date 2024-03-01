package GenKode::BlastParse;
use strict;
use warnings;

use File::Spec;
use File::Path;
use File::Temp;
use Data::Dumper;
use LWP::Simple;
use Storable qw/dclone nstore/;
use XML::Simple;

#for remote blast
use URI::Escape;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);

# ----------- module for parsing ncbi's blast output ------------

# 200207: xml to txt
# 181219: tested for local blast against blast-db . 1 query per xml

sub new { # feed in one valid blast xml output file
# xml hash is stored in $self->{xml} <<< A ref
	my ($class, $file)=@_;
	my $self= bless { searches =>[] },$class;
	eval {
		my ($xx)=_blast2href($file);
		$self->{searches}=dclone($xx);
		# $self->{query}=dclone ($qq);
		# $self->{meta}=dclone ($mm);
	};
	return $self;
}
sub _blast2href { #only keep hits and param info. works for both xml=5 or 14
	my ($file)=@_;
	if (!$file or !-e $file) {
		return ({},[]);
	}

	# my $query;
	# if (open my $fh, $file) { #for query, quicker to just open the file
		# while (<$fh>) {
			# chomp;
			# if (/<query-(len|id|title)>(.+?)<\/query-.+?>/) {
				# $query->{$1}=$2;
			# }
		# }
		# close ($fh);
	# }

	my $report;
	my $xr0=_xml2hash($file,['BlastOutput2','HitDescr','Hit','Hsp'], 1);
	my $xa=[];
	#ONLY xml contain ONE iteration is processed currently (19apr18)
	if ($xr0->{BlastOutput2}) { # GB server has this
		foreach my $search (@{$xr0->{BlastOutput2}}) {
			my $mref=$search->{report}{Report}{results}{Results};
			$report=$search->{report}{Report};
			if ($mref->{iterations}{Iteration}) {
				push @$xa, $mref->{iterations}{Iteration}{search}{Search};
			} else {
				push @$xa, $mref->{search}{Search};
			}
		}
	}
	else { # blast to seq has bl2seq
		$report=$xr0->{report}{Report};
		foreach my $m (qw/ bl2seq  search /) {
			my $v = $xr0->{report}{Report}{results}{Results}{$m}{Search};
			push @$xa, $v if defined $v;
		}
	}

	# my $meta;
	# foreach my $k (qw/program version reference /) { #supposed to be string
		# $meta->{$k}=$report->{$k}
	# }
	# if ($report->{params}{Parameters}) {
		# $meta->{params}=dclone $report->{params}{Parameters};
	# }
	# return ($query,$xa,$meta);
	return ($xa);
}
sub _xml2hash {
	my ($file,$force, $keyattroff)=@_; #file, A ref
	my $xml=XML::Simple->new();
	my $ref;
	eval {
		if ($keyattroff) {
			$ref=$xml->XMLin($file, ForceArray => $force, KeyAttr => [ ]);
		} else {
			$ref=$xml->XMLin($file, ForceArray => $force);
		}
	};
	if ($@) {
		return {};
	} else {
		return $ref;
	}
}

sub _chk_idx { #verify if input is a number. MUST be positive. >0
	my $n=shift;
	if ($n and $n=~/^\d+$/) {
		return 1;
	} else {
		return 0;
	}
}

sub nomatch { #check whether "no hits found"
	my ($self,$sid)=@_;
	$sid=1 if !$sid;
	my $sr=$self->get_search($sid);
	if ($sr->{message} and $sr->{message}=~/no hits/i) {
		return 1;
	} else {
		return 0;
	}
}

sub all_searches { # return number of searches
	my $self=shift;
	return scalar @{$self->{searches}}||0;
}
sub get_search { #by index
	my ($self,$sid)=@_; #if idx given, return the specified H ref, else return an A ref of all idx
	if (_chk_idx($sid)) {
		$sid--; #idx is in human order
		return $self->{searches}[$sid]||{};
	}
	return {};
}
sub all_hits { #return number of hits belong to the specified search-index
	my ($self,$sid)=@_;
	my $xs=$self->get_search($sid);
	if (defined $xs->{hits}{Hit}) {
		return scalar @{$xs->{hits}{Hit}}||0;
	}
	return 0;
}
sub get_hit { #by giving search-index and hit-index
	my ($self,$sid,$hid)=@_;
	my $xs=$self->get_search($sid);
	if ($xs and _chk_idx($hid)) {
		$hid--;
		return $xs->{hits}{Hit}[$hid]||{};
	}
	return {};
}
sub all_hsps {
	my ($self,$sid,$hid)=@_;
	if (my $hit=$self->get_hit($sid,$hid)) { #if hit exists, hsp must contain at least 1 entry
		return scalar @{$hit->{hsps}{Hsp}}||0;
	}
	return 0;
}
sub get_hsp {
	my ($self,$sid,$hid,$pid)=@_;
	my $hit=$self->get_hit($sid,$hid);
	if ($hit and _chk_idx($pid)) {
		$pid--;
		return $hit->{hsps}{Hsp}[$pid]||{};
	}
	return {};
}

#----only work on one single hsp.----
#$hsp should be from $hsp=$self->get_hsp($sid,$hid,$pid)
sub identity { # for aligned region only
	my ($self,$hsp)=@_;
	return ($hsp->{identity}/$hsp->{'align-len'})||-1;
}
sub coverage { # (query_from - query_to + 1) / query_len
	my ($self, $hsp, $sid)=@_;
	# return ($hsp->{'align-len'}/$self->query('len', $sid));
	return ((abs($hsp->{'query-from'} - $hsp->{'query-to'}) +1 )/$self->query('len', $sid));
}

sub query { #known:  query-masking -len -id -title
	my ($self,$var, $sid)=@_; #return string if "var" specified . return H ref if no var
	$sid=1 if !$sid;
	my $search=$self->get_search($sid);
	if ($var) {
		if ($var=~/len(gth)?/i) {
			return $search->{'query-len'};
		} else {
			$var='query-'.(lc $var);
			return $search->{$var}||'';
		}
	} else {
		return '';
	}
}
=pod
sub meta {
	my ($self,$var,$isparam)=@_; #return string if "var" specified . program|version|reference || meta: expect|sc-match|sc-mismatch|gap-open|gap-extend|filter
	if ($isparam) {
		return $self->{meta}{params}{$var}||'';
	}
	elsif ($var) {
		$var=lc $var;
		return $self->{meta}{$var}||'';
	}
	else {
		return $self->{meta};
	}
}
=cut

sub hit { # title id accession taxid sciname?
	my ($self,$sid,$hid,$var)=@_;
	my $hit=$self->get_hit($sid,$hid);
	$var='' if !$var;
	if ($var) {
		$var=lc $var;
		if ($hit) {
			if ($var=~/len(gth)?/) {
				return $hit->{len};
			} else { #for shortened names
				if ($var=~/^acc/) {
					$var='accession';
					my $id=$hit->{description}{HitDescr}[0]{'id'};
					my $acc1=$hit->{description}{HitDescr}[0]{'accession'};
					my $acc=$acc1;
					if ($id=~/$acc1(\.\d+)/i) {
						$acc.=$1;
					}
					return $acc;
				}
				elsif ($var=~/^tax(on)?/) { $var='taxid'; }
				return $hit->{description}{HitDescr}[0]{$var}||'';
			}
		}
	} else {
		return $hit->{description}{HitDescr}[0];
	}
}
sub hsp { # num bit-score score evalue identity query-from query-to query-strand hit-from hit-to hit-strand align-len gaps qseq hseq midline
	my ($self,$sid,$hid,$pid,$var)=@_;
	$var='' if !$var;
	$var=lc $var;
	if (my $hsp=$self->get_hsp($sid,$hid,$pid)) {
		return defined($hsp->{$var})?$hsp->{$var}:'';
	}
	return '';
}

#---- extract blast-based alignment, mostly non-OO -----
sub extract_aln_shell {
	my ($self, $fh2, $query, $xmls, $filters)=@_;
	$fh2=*STDOUT if !$fh2;
	$query=~s/\W//g;
	my $qlen=length $query;
	my $xalns=$self->extract_aln($xmls,$filters);
	my $colmax=$self->calc_col_len($xalns, $qlen);
	$self->print_aln_query($query,$colmax,$fh2);
	foreach my $kk (keys %$xmls) {
		printf $fh2 ">file__%s\nxxxxx\n", $kk;
		$xmls->{$kk}->print_aln_hsp_all($xalns->{$kk},$colmax,$fh2);
	}
	1;
}
sub extract_aln { # input multiple $xr objs, return hash for all aln-info by position. return A ref
	my ($self,$xmls, $filters)=@_; # filter key same as xmls key
	my $xalns;
	foreach my $fn (keys %$xmls) {
	# die Dumper $xmls->{$fn},$fn;
		my $aln=$xmls->{$fn}->extract_aln_single($filters->{$fn});
		$xalns->{$fn}=dclone $aln;
	}
	return $xalns;
}
sub extract_aln_single { #return alnset data. if one xml contains multiple searches, only take the first. all xml s are supposed to be from the same query
	my ($xr,$filter)=@_; #filter given as H ref { hit=>{hsp=>1, hsp2=>1, ...}, hit2=>{} , ... }
	my $alnset;
	# die Dumper [keys %$filter];
	my $sid=1; #only process 1st search
	foreach my $hid (1..$xr->all_hits($sid)) {
		#-------- get alignment
		if ($filter and !defined $filter->{$hid}) {
			next;
		}
		foreach my $tid (1..$xr->all_hsps($sid,$hid)) { #loop per hsp
			if ($filter and !defined $filter->{$hid}{$tid}) {
				next;
			}
			# printf "hit=%s, hsp=%s\n", $hid, $tid;
			my $harr=$xr->extract_aln_by_hsp($sid,$hid,$tid);
			$alnset->{$hid}{$tid}=dclone $harr;
		}
	}
	return $alnset;
}
sub extract_aln_by_hsp { #this allows to process any selected hsp in any search if multiple are included. but need to make the loop outside of this package
	my ($xr,$sid,$hid,$tid)=@_;
	my $hsp=$xr->get_hsp($sid,$hid,$tid);
	my $currq=$hsp->{'query-from'};
	my $harr;
	for (my $j=0; $j<length $hsp->{'midline'}; $j++) {
		my $cq=substr $hsp->{qseq}, $j, 1;
		my $ch=substr $hsp->{hseq}, $j, 1;
		if ($cq=~/\w|\*/) { #cq isn't a gap
			$harr->{$currq}=$ch;
			$currq++;
		} else { #connect hit AA to previous
			if (!$harr->{$currq-1}) {
				$harr->{$currq-1}='';
			}
			$harr->{$currq-1}.=$ch;
		}
	}
	return $harr;
}
sub calc_col_len { # input &extract_aln_multi processed aln A ref data, calculate the max col-span for each character
	my ($self, $xalns,$qlen)=@_;
	my $colmax;
	foreach my $fn (keys %$xalns) {
		my $aln_hit=$xalns->{$fn};
		foreach my $hid (keys %{$aln_hit}) {
			foreach my $tid (keys %{$aln_hit->{$hid}}) {
				my $aln1=$aln_hit->{$hid}{$tid};
				for my $pos (1..$qlen) {
					if ($aln1->{$pos}) {
						my $plen=length $aln1->{$pos};
						if (!$colmax->{$pos} or $plen>$colmax->{$pos}) {
							$colmax->{$pos}=$plen;
						}
					} else {
						if (!$colmax->{$pos}) {
							$colmax->{$pos}=1;
						}
					}
				}
			}
		}
	}
	$colmax->{maxlen}=$qlen;
	return $colmax;
}
sub print_aln_query { # print query for extracted aln
	my ($self,$query,$colmax,$fh2)=@_; #must be full length query as used for blast
	printf $fh2 ">query\n";
	for my $pos (1..$colmax->{maxlen}) { #foreach position
		if (!$colmax->{$pos}) {
			$colmax->{$pos}=1;
		}
		my $collen=$colmax->{$pos};
		print $fh2 (substr $query, ($pos-1), 1);
		if ($collen>1) {
			print $fh2 ('-' x ($collen-1));
		}
	}
	print $fh2 "\n";
}
sub print_aln_hsp_all {
	my ($self,$aln_hit,$colmax,$fh2)=@_;
	# die Dumper [keys %$aln_hit];
	foreach my $hid (sort {$a<=>$b} keys %{$aln_hit}) {
	# die Dumper $aln_hit->{$hid};
		foreach my $tid (sort {$a<=>$b} keys %{$aln_hit->{$hid}}) {
			my $aln=$aln_hit->{$hid}{$tid};
			# die Dumper $aln;
			# printf "    hit=%s, hsp=%s\n", $hid, $tid;
			print_aln_hsp_one($aln, $colmax, (sprintf 'hit_%i__hsp_%i', $hid,$tid), $fh2);
		}
	}
}
sub print_aln_hsp_one { #to print one seq. flexible for only need specific seqs
	my ($aln, $colmax, $sname, $fh2)=@_; #aln is sub key contains ONE  hit/hsp info.. maxlen can use length of query
	$fh2=*STDOUT if !$fh2;
	$sname=sprintf 'hsp_%s', time if !$sname;
	printf $fh2 ">%s\n", $sname;
	for my $pos (1..$colmax->{maxlen}) { #foreach position
		if (!$aln->{$pos}) { #nothing at that position
			print $fh2 ("-" x $colmax->{$pos});
		} else {
			print $fh2 $aln->{$pos};
			if ((length $aln->{$pos}) < $colmax->{$pos}) {
				print $fh2 ('-' x ($colmax->{$pos} - (length $aln->{$pos})));
			}
		}
	}
	print $fh2 "\n";
}

# --------- convert xml2 format to plain text format ------------
sub print_aln_block {
	my ($xr,$sid,$hid,$pid,$fh,$num)=@_;
	$num=60 if !$num;
	$fh=*STDOUT if !$fh;
	my $hsp=$xr->get_hsp($sid,$hid,$pid);
	my $qs=$hsp->{qseq};
	my $hs=$hsp->{hseq};
	my $ms=$hsp->{midline};
	my $lenm=0;
	foreach my $k (qw/query-from  hit-from /) {
		if ($lenm<(length $hsp->{$k})) {
			$lenm=length $hsp->{$k};
		}
	}
	if ($num eq '-1') {
		printf $fh "Query %*d  %s  %s\n", $lenm, $hsp->{'query-from'}, $qs, $hsp->{'query-to'};
		printf $fh "      %s  %s\n", (' ' x $lenm), $ms;
		printf $fh "Sbjct %*d  %s  %s\n\n\n", $lenm, $hsp->{'hit-from'}, $hs, $hsp->{'hit-to'};
	}
	else {
		my $cend=0;
		my $qbegin=$hsp->{'query-from'};
		my $hbegin=$hsp->{'hit-from'};
		my $plus=$hsp->{'hit-strand'} eq 'Plus'?1:-1;
		my $dist=$hsp->{'query-to'}-$hsp->{'query-from'}+1;

		while ($cend<$dist) {
			if (($cend+$num)>$dist) {
				$num=$dist-$cend;
			}
			my $q1=substr $qs, $cend, $num;
			my $h1=substr $hs, $cend, $num;
			my $m1=substr $ms, $cend, $num;
			printf $fh "Query %*d  %s  %s\n", $lenm, $qbegin, $q1, ($qbegin+$num-1);
			printf $fh "      %s  %s\n", (' ' x $lenm), $m1;
			printf $fh "Sbjct %*d  %s  %s\n\n", $lenm, $hbegin, $h1, ($hbegin+($num-1)*$plus);
			$qbegin+=$num;
			$hbegin+=$num;
			$cend+=$num;
		}
		printf $fh "\n";
	}
}
sub print_out_head {
	my ($xr,$fh,$sid)=@_;
	$fh=*STDOUT if !$fh;
	printf $fh "%s\n\n\n", $xr->meta('version');
	printf $fh "%s\n\n\n", $xr->meta('reference');
	printf $fh "Query= %s\n\n", $xr->query('title');
	printf $fh "Length= %s\n\n", $xr->query('len');
	printf $fh "Sequences producing significant alignments (tab-delimited)\n\n";
	printf $fh "Title\tScore(bits)\tE-value\n";
	#summary chart
	foreach my $hid (1..$xr->all_hits($sid)) { #suppose only has 1 search
		printf $fh "%s...\t%.1f\t%.0e\n",
			substr ($xr->hit($sid,$hid,'title'), 0, 65),
			$xr->hsp($sid,$hid,1,'bit-score'),
			$xr->hsp($sid,$hid,1,'evalue');
	}

	print $fh "\n\n\n";
}
sub xml2txt { #directly print to a fh, (near-)standard blast xml format
	my ($xr, $fh, $param)=@_;
	$fh=*STDOUT if !$fh;
	my $num=$param->{num}||60;
	if (!$param->{simple}) {
		$xr->print_out_head($fh,$param->{sid}||1);
	}
	#matched alignment
	foreach my $sid (1..$param->{search}||$xr->all_searches) {
		foreach my $hid (1..$param->{hit}||$xr->all_hits($sid)) {
			printf $fh ">%s\n",$xr->hit($sid,$hid,'title');
			if (!$param->{simple}) {
				printf $fh "Length=%s\n",$xr->hit($sid,$hid,'len');
			}
			foreach my $pid (1..$param->{hsp}||$xr->all_hsps($sid,$hid)) {
				print $fh "\n";
				my $hsp=$xr->get_hsp($sid,$hid,$pid);
				if (!$param->{simple}) {
					printf $fh "  Score = %s bits (%s),  Expect = %s\n", $hsp->{'bit-score'}, $hsp->{score}, $hsp->{evalue};
					printf $fh "  Identities = %s/%s (%.0f%%), Gaps = %s/%s (%.0f%%)\n",
						$hsp->{identity}, $hsp->{'align-len'}, ($hsp->{'identity'}/$hsp->{'bit-score'}),
						$hsp->{'gaps'}, $hsp->{'align-len'}, ($hsp->{'gaps'}/$hsp->{'align-len'});
					printf $fh "  Strand=Plus/%s\n\n", $hsp->{'hit-strand'};
				}
				$xr->print_aln_block($sid,$hid,$pid,$fh,$num);
			}
		}
	}
}

1;

