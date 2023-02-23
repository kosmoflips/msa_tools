package GenKode::Fold;
use strict;
use warnings;

# specialised to only process folding formats: vna(fst), ct, mfold constraints

use Storable qw/dclone nstore retrieve/;
use File::Spec;
use File::Path;
use Data::Dumper;

use base 'Exporter';
our @EXPORT = qw/
mkfold
/;

our @EXPORT_OK = qw/

/;

our %EXPORT_TAGS = (
# translation=>[qw//]
);


# https://perldoc.perl.org/5.8.8/Class/ISA.html
# use ISA to access some MSA methods

use Class::ISA; #install if not found!

use GenKode::MSA;
@GenKode::Fold::ISA=qw(GenKode::MSA);
# print Class::ISA::super_path('GenKode::Fold');
# exit;



sub new {
	my ($class)=shift @_;
	my $self= bless {},$class;
	$self->{env}{_lastid}=0;
	$self->{env}{_maxlen}=0;
	$self->{env}{_maxlen_name}=0;
	return $self;
}
sub import { #import a file to the Obj. only work for fst/vna/ct. must have a seq and a fold. !! DOESN'T input mfold constraints, but can output folding info in Mfold constraints format
	my ($self, $filepath) = @_;
	if (!$filepath or !-e $filepath or -B $filepath) {
		return 0; # "can't properly open file or file doesn't exist";
	}
	my $format= $self->_guess_fold($filepath);
	open (my $fh, $filepath);
	if ($format eq 'ct') {
		my @c=File::Spec->splitpath ($filepath);
		my $fname=pop @c||'';
		$fname=~s/\.ct$//i;
		$self->read_fh_ct($fh,$fname);
	}
	elsif ($format eq 'vna') {
		$self->read_fh_vna($fh);
	}
	elsif ($format eq 'fst') {
		$self->read_fh_fst($fh); #from MSA
		$self->parse_fold;
	}
	else {
		return 0; # "unknown file format";
	}
	$self->_maxlen_batch;
	$self->_maxnamelen_batch;
	return 1;
}

#simple convert between formats. if want more options (get single seq per file, get masked region, use MSA instead)
sub _eval_outformat {
	my $format=shift;
	if ($format=~/fst|fasta|fsa/i) { return 'fst'; }
	elsif ($format=~/mfold/i) { return 'mfold'; }
	elsif ($format=~/ct/i) { return 'ct'; }
	elsif ($format=~/vna|vienna/i) { return 'vna'; }
	return '';
}
sub export { #fst , vna , ct, mfold
	my ($self, $infile, $format, $fixfname, $voff)=@_; #$voff is only for write-to-fst
	$format=_eval_outformat($format);
	if ($format eq 'fst') {
		my $ofile=$fixfname?$infile:(sprintf '%s%s.fst', $infile,($voff?'_voff':''));
		open (my $fh2, ">", $ofile) or return 0;
		$self->write_fh_fst($fh2,$voff);
		return $ofile;
	}
	elsif ($format =~/vna|ct/) {
		my $odir;
		if ($self->lastfold>1) {
			$odir=$fixfname?$infile:(sprintf '%s_Fold_%s', $infile,$format);
			mkpath $odir;
			if (!-d $odir) {
				return 'can\'t make folder for output';
			}
		}
		foreach my $id (1..$self->lastid) {
			next if !$self->has_vna($id); #this sub also chk if it exists
			my $sname=GenKode::MSA::safename($self->getname($id));
			my $ofile;
			if ($odir) {
				$ofile=File::Spec->catfile($odir, $sname.'.'.$format);
			} else {
				$ofile=$fixfname?$infile:$infile.'.'.$format;
			}
			open (my $fh2, ">", $ofile) or next;
			if ($format eq 'vna') {
				$self->write_fh_vna($fh2,$id);
			} else {
				$self->write_fh_ct($fh2,$id);
			}
			if ($self->lastfold==1) { #only 1 folding vna, can return
				return $ofile;
			}
		}
		return $odir;
	}
	elsif ($format eq 'mfold') {
		my $ofile=$fixfname?$infile:(sprintf '%s_mfold.txt', $infile);
		open (my $fh2, ">", $ofile) or return 0;
		$self->write_fh_mfold($fh2);
		return $ofile;
	}
	return 'nothing exported. .';
}
sub export_error {
	my ($self, $infile, $fixfname)=@_;
	my $errfile=$fixfname?$infile:$infile.'_errors.txt';
	open (my $fh3, ">", $errfile);
	my $has_err;
	foreach my $id (1..$self->lastid) {
		next if !$self->has_vna($id);
		my $fold=mkfold($self->getseq($id),$self->getvna($id));
		if ($fold->[0]) {
			printf $fh3 "--[%s] %s--\n", $id,$self->getname($id);
			$self->write_fh_error($fh3, $fold->[0]);
			print $fh3 "//\n";
			$has_err=1;
		}
	}
	if ($has_err) {
		return $errfile;
	} else {
		close ($fh3);
		unlink $errfile;
		return 0;
	}
}

sub fold2seq { #non-OO
	my ($fold)=@_; #A ref
	my $s='';
	my $v='';
	shift @$fold;
	foreach my $r (@$fold) {
		$s.=$r->[0];
		$v.=$r->[2]||'.';
	}
	return ($s,$v);
}
sub map_fold { #match a given fold vna to seq by id. seq must be exact match
#NOTE: didn't verify how this sub works if seqs do not exactly match
	my ($self, $id, $fold)=@_; #fold is A ref made by mkfold
	my $ts=$self->getseq($id); #template seq
	return '' if !$fold;

	my ($fseq,$fvna)=fold2seq($fold);
	my $ts0=$ts; $ts0=~s/\W//g; $ts0=~s/t/u/ig;
	my $fseq0=$fseq; $fseq0=~s/t/u/ig;
	if (!$fseq or !$fseq or $ts0!~/$fseq0/i) {
		return '';
	}
	$ts0=~/$fseq/i;
	my ($pre,$post)=($`,$');
	my $v='';
	my $j=0;
	my $r=0; #real letter
	for (my $i=0; $i<length $ts; $i++) {
		my $ct=substr $ts, $i, 1;
		if ($ct=~/[a-z]/i) {
			$r++;
			next if $r<=length $pre;
			last if $j>length $fseq;
			my $cf=substr $fvna, $j, 1;
			$v.=$cf;
			$j++;
		} else {
			$v.='.';
		}
	}
	$self->{$id}{vna}=$v;
	return $v;
}

# ------- read files ---------
sub _guess_fold { #guess input format , among vna/fst/ct. vna format should only contain 1 seq
	my ($self, $path)=@_;
	open (my $fh, $path);
	return 0 if ref $fh !~/glob/i;

	my $line=0;
	my $vna=0; #count of vna folding line ...((...))...
	my $fst=0; #count of fasta title /^>/
	my $is_ct=0;
	my $fsttitle=0;
	my $seqline=0;
	my $vnaline=0;
	while (<$fh>) {
		chomp;
		$line++;
		if ($line>1 and
			/^ \s* \d+  \s+    \w  \s+   \d+ \s+   \d+ \s+   \d+ \s+   \d+  \s*$/x
		) { #not the 1st line. can guess ct
			#8	g	7	9	210	8
			$is_ct=1;
			last;
		}
		#guess vna: line1=title,line2=seq,line3=vna,line4=title,line5=seq....
		# if ($line%3==0 and /^[^>]/ and /^[^a-z]+$/i ) { $vna++; }
		if (/^>/ and $line==1) { $fsttitle++;}
		elsif ($line==2 and /[a-z]/i) { $seqline++; }
		elsif ($line==3 and /^[^a-z]+$/i) { $vnaline++; }
		# elsif (/^[^a-z]+$/i) {
			# $vnaline++;
		# }
		# elsif (/[a-z]/i) {
			# $seqline++;
		# }
	}
	if ($seqline==$vnaline and $line==3) {
		return 'vna';
	}
	# if ($fsttitle and $seqline==$vnaline and $seqline==$fsttitle) {
		# return 'vna';
	# }
	elsif ($fsttitle>0) {
		return 'fst';
	}
	elsif ($is_ct) {
		return 'ct';
	}
	else {
		return 0; #unknown
	}
}
sub read_fh_ct {#ct, supposed to contain only one sequence. REMEMBER TO FEED IN A NAME.
# pos NT prev next pair pos
	my ($self,$fh,$name)=@_; #input is a FH
	my ($seq,$vna)=('','');
	while (<$fh>) {
		chomp;
		if (!$name and /^\s*(\d+)\s+/) { #use whatever on 1st line as the name. ^\s*\d+\s+(.+)
			$name=$';
			$name='seq_'.time() if !$name;
		}
		next if !/^\s*\d+\s+[atcgunx]/i; # skip non standard lines. shouldn't happen in a well formatted ct file
		my $line=$_;
		$line=~s/^\s+//;
		my @col=split /\s+/, $line;
		next if @col!=6; #shouldn't happen
		$seq .=$col[1];
		if ($col[0]<$col[4]) {
			$vna .= '(';
		} elsif ($col[4]==0) {
			$vna .= '.';
		} else {
			$vna .=')';
		}
	}
	$self->lastid(1);
	$self->lastfold(1);
	$self->{$self->lastid}=dclone {seq=>$seq, vna=>$vna, name=>$name};
	1;
}
sub read_fh_vna { #vna format: line 1,2,3: >name, seq, vna ... otherwise won't work. does NOT include any special seqs (e.g. annotation line)
	my ($self,$fh)=@_;
	my $nameafter;
	my $idxbegin=$self->lastid+1;
	while (<$fh>) {
		chomp;
		my $line=$_;
		$line=~s/\s//g;
		next if /^\s*$/;
		if ($line=~/^>/) {
			$self->lastid(1);
			$self->{$self->lastid}{name}=$';
			$nameafter=1;
		} elsif ($nameafter) {
			$nameafter=0;
			$self->{$self->lastid}{seq}=$line;
		} else {
			$self->{$self->lastid}{vna}=$line;
		}
	}
	$self->lastfold(1);
	1;
}
# sub read_fh_fst . from MSA
sub parse_fold { #fst is processed from MSA method, need to combine vna info // tricky since annotation lines can also have {}[]()
	my ($self)=@_;
	my $seqok;
	$self->{env}{_lastfold}=0;
	foreach my $id (1..$self->lastid) {
		next if !$self->contains($id);
		my $stype=_eval_seq_type($self->getseq($id));
		if ($stype==1) {
			$seqok=$id;
		}
		elsif ($seqok and $stype==2 and ($seqok+1)==$id) { #this folding must be right after seq
			$self->{$seqok}{vna}=$self->getseq($id);
			$self->lastfold(1);
			$self->drop($id);
			$seqok=0;
		}
		# print $self->getseq($id),$stype;<>;
	}
	1;
}
sub lastfold {
	my ($self,$raise)=@_;
	$self->{env}{_lastfold}=0 if !$self->{env}{_lastfold};
	if ($raise) {
		$self->{env}{_lastfold}++;
	}
	return $self->{env}{_lastfold};
}
sub _eval_seq_type { # non-OO
	my $seq=shift;
	if (!$seq) {
		return 0;
	}

	my $pair=0;
	my $digit=0;
	my $letter=0;
	for my $i0 (1..length $seq) {
		my $i=$i0-1;
		my $c=substr $seq, $i, 1;
		if ($c=~/[a-z\*]/i) { $letter++; }
		elsif ($c=~/\d/) { $digit++; }
		elsif ($c=~/[  \[\]\(\)\{\}  ]/x) { $pair++; }
	}

	if ($pair<$letter) { #more letters and star(stop codon) only.#likely a real nt/aa seq
		return 1; 
	}
	elsif ($pair>$letter) {#contain pairing symbols
		return 2; #likely a folding vna seq
	}
	#elsif .... for AA folding annotation?
	return 0;
}

sub write_fh_fst {
	my ($self,$fh,$voff,$id)=@_;
	foreach my $i (1..$self->lastid) {
		if ($id) {
			next if $i!=$id;
		}
		next if !$self->contains($i);
		$self->write_fh_fst_1seq($fh, $self->{$i}{name}, $self->{$i}{seq});
		if (!$voff and $self->has_vna($i) ) {
			$self->write_fh_fst_1seq($fh, $self->{$i}{name}.'_fold', $self->{$i}{vna});
		}
	}
	1;
}
sub _pureseq { #get pure sequence without gap. access {ct}
	my ($self, $i, $usedot)=@_;
	my $s='';
	my $v='';
	if ($self->has_vna($i)) {
		foreach (my $j=0; $j< length $self->{$i}{seq}; $j++) {
			my $c1=substr $self->{$i}{seq}, $j, 1;
			my $c2=substr $self->{$i}{vna}, $j, 1;
			if (_is_gap($c1) + _is_gap($c2) < 2) {
				$s.=$c1;
				if ($c2=~/[a-z]/i or _is_gap($c2)) {
					$v.=$usedot?'.':'-';
				} else {
					$v.=$c2;
				}
			}
		}
	}
	return ($s,$v);
}
sub _is_gap {
	my $string=shift;
	if ($string=~/^[\s\_\-\.~\,]*$/) {
		return 1;
	} else {
		return 0;
	}
}
sub write_fh_vna {#write ONE seq only
	my ($self,$fh,$i)=@_;
	my ($s,$v)=$self->_pureseq($i,1);
	printf $fh ">%s\n%s\n%s\n", $self->getname($i), $s,$v;
	1;
}
sub write_fh_fold { #fold A ref is ready
	my ($self, $fh, $fold, $name)=@_;
	printf $fh "%i\t%s\n", ( (scalar @{$fold}) -1 ), $name||'sequence_fold';
	for (my $idx=1; $idx < scalar @{$fold} ; $idx++) {
		my $ref=$fold->[$idx];
		printf $fh "%s\t%s\t%s\t%s\t%s\t%s\n",$idx,$ref->[0],$idx-1,$idx+1,($ref->[1]||0),$idx;
	}
	1;
}
sub write_fh_ct {
	my ($self,$fh,$i)=@_;
	my $fold=mkfold($self->{$i}{seq},$self->{$i}{vna});
	$self->write_fh_fold($fh, $fold, $self->getname($i));
	1;
}
sub write_fh_mfold { # note, if folding vna contains pseudoknot, Mfold won't work
	my ($self,$fh)=@_;
	print $fh "*Note! if your folding contains pseudoknots, generated constraints won't work on Mfold (you'll be informed by Mfold)!\n\n";
	foreach my $id (1..$self->lastid) {
		next if !$self->has_vna($id);
		$self->write_fh_mfold_1seq($fh,$id);
	}
	1;
}
sub write_fh_mfold_constrain { #constrain block only
	my ($self,$fh,$fold,$content)=@_; #content: 0: write both/1: write force pair only "F i j k" /2: write single strand only "P i j k"
	my $pair; my $ss;
	my $last_p=0;
	my $last_q=0;
	my $ss_p=0;
	my $ss_q=0;
	for (my $p=1; $p< scalar @{$fold} ; $p++) {
		my $q=$fold->[$p][1] || 0;
		my $step=$pair->{$last_p}[1]||0;
		if ($q) {
			if ($p<$q) { #processing pairs. p is paired to q
				if ( ( ($last_p+ $step ) == $p ) and 
				( ($last_q- $step ) == $q ) ) { #last p, q are next to current p,  q
					if ($pair->{$last_p}[1]) {
						$pair->{$last_p}[1]++;
					} else {
						$pair->{$last_p}[1]=1;
					}
				} else {
					$last_p=$p;
					$last_q=$q;
					$pair->{$last_p}=[$last_q, 1];
				}
			} else {
				$last_p=0;
				$last_q=0;
			}
		} else {
			if (!$ss_p) {
				$ss_p=$p;
				$ss_q=$p;
			}
			elsif ($ss_q+1 == $p) {
				$ss_q=$p;
			}
			else {
				$ss->{$ss_p}=$ss_q-$ss_p+1;
				$ss_p=$p;
				$ss_q=$p;
			}
		}
		#process trailing ss region
		if ($ss_p) {
			$ss->{$ss_p}=scalar @{$fold} - $ss_p; #scalar contains elem[0], so no need to +1
		}
	}

	$content=0 if !$content;
	if ($content==0 or $content==1) { #pair
		foreach my $ff (sort {$a<=>$b} keys %$pair) {
			next if $ff==0;
			#force consecutive base pairs i.j,i+1.j-1, ...,i+k-1.j-k+1
			if ($pair->{$ff}[1] > 1) { #to avoid single basepair
				printf $fh "F %i %i %i\n", $ff, $pair->{$ff}[0], $pair->{$ff}[1];
			}
		}
	}
	if ($content==0 or $content==2) { #ss
		foreach my $pp (sort {$a<=>$b} keys %$ss) {
			next if $pp==0;
			printf $fh "P %i 0 %i\n", $pp, $ss->{$pp};
		}
	}
	1;
}
sub write_fh_mfold_1seq {
	my ($self,$fh,$i)=@_;
	my ($s,$v)=$self->_pureseq($i);
	$self->write_fh_vna($fh, $i);
	print $fh "\n# constraints\n";
	my $fold=mkfold($s,$v);
	if ($fold->[0]) {
		printf $fh "!! there are errors in given file. fix them first\n";
		$self->write_fh_error($fh,$fold->[0]);
	}
	else {	# process F/P info
		$self->write_fh_mfold_constrain($fh,$fold);
	}
	printf $fh "//\n";
}
sub write_fh_error { #for one seq
	my ($self, $fherr, $errinfo)=@_;
	if ($errinfo) { # A ref
		foreach my $type (qw/left-open right-open mismatch N-pair/) { ##assigned in mkfold
			next if !$errinfo->{$type};
			printf $fherr "==%s==\n", $type;
			foreach my $el (@{$errinfo->{$type}}) {
				printf $fherr "%s\n", $el;
			}
		}
		return 1;
	} else {
		return 0;
	}
}
sub mkfold_id { #for OO use
	my ($self,$id)=@_;
	return mkfold($self->getseq($id),$self->getvna($id));
}
sub mkfold { # valid pairs: ( ) [ ] { } . do not mix them, ( ] isn't a valid pair. [0] is for error
	my ($seq,$vna)=@_;
	my $symbol={
		'(' => ')',
		'{' => '}',
		'[' => ']',
	};
	my $half; #positions of opening "(  [  {" labelled
	my $fold;
	my $ferr;
	my $char=0;
	my $i=0; #human order along the sequence

	for my $char (0.. ((length $seq)-1) ) {
		my $nt=substr $seq,$char,1;
		if ($nt=~/[a-z]/i) { # must be a valid NT. so excluding "-" "." etc..
			$i++;
			$fold->[$i][0]=$nt;
			my $v='-';
			if ($char<length $vna) {
				$v=substr $vna,$char,1;
			}
			my $is_paired=0;
			foreach my $ss (keys %$symbol) {
				if ($v eq $ss) { #[ ( {
					$is_paired=1;
					push @{$half->{$ss}}, $i;
				}
				elsif ($v eq $symbol->{$ss}) { #] ) }
					$is_paired=1;
					if (!$half->{$ss} or
						( (ref $half->{$ss} eq 'ARRAY') and (scalar @{$half->{$ss}} == 0) ) ) { #missing opening half
						push @{$ferr->{'left-open'}}, (sprintf "? -> %s%i",$fold->[$i][0],$i);
					} else {
						my $pair=pop @{$half->{$ss}};
						$fold->[$i][1]=$pair;
						$fold->[$pair][1]=$i;
						my $e=eval_pair($fold->[$pair][0],$fold->[$i][0]);
						if ($e!=1 and $e!=2) {
							my $mark=($e==-1?'N-pair':'mismatch');
							push @{$ferr->{$mark}}, (sprintf "%s%s <?> %s%s", $fold->[$pair][0],$pair, $fold->[$i][0],$i);
						}
					}
				}
			}
			if (!$is_paired) {
				$fold->[$i][1]=0;
			} else {
				$fold->[$i][2]=$v;
			}
		}
	}
	foreach my $ss (keys %{$symbol}) {
		if ($half->{$ss}) { #perfectly paired sequence should have no remaining @half elements
			foreach my $pos (@{$half->{$ss}}) {
				push @{$ferr->{'right-open'}}, (sprintf "%s%i <- ?",$fold->[$pos][0],$pos);
			}
		}
	}
	$fold->[0]=$ferr?(dclone $ferr):undef;
	return $fold;
}
sub eval_pair { #non-OO // for ct. return: 1:good WC pairing. 2:wobble pairing. -1:ambiguous N to N. 0:bad
	my ($nt,$nt2)=@_;
	$nt=lc $nt;
	$nt2=lc $nt2;
	if ($nt eq 'u') { $nt='t'; }
	if ($nt2 eq 'u') { $nt2='t'; }
	my $eval=0;
	if (
		($nt eq 'a' and $nt2 eq 't') or
		($nt eq 't' and $nt2 eq 'a') or
		($nt eq 'c' and $nt2 eq 'g') or
		($nt eq 'g' and $nt2 eq 'c')
	) {
		$eval=1;
	}
	elsif (
		($nt eq 't' and $nt2 eq 'g') or
		($nt eq 'g' and $nt2 eq 't')
	) {
		$eval=2;
	}
	elsif ($nt=~/n/i or $nt2=~/n/i) {
		$eval=-1;
	}
	return $eval;
}

1;