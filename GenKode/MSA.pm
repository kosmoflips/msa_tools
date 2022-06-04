package GenKode::MSA;
use strict;
use warnings;

use Storable qw/dclone/;
use File::Spec;
use File::Path;
use Data::Dumper;


use GenKode::CentralDogma;

### manipulation MSA related things, format converts, translation, etc.
### by kiyo @ http://www.pocchong.de
### created: 2013-10-05
# ----- 181211 ----
# with some major revision
# each seq is treated independently. do NOT consider/combine/process structure/fold/ct, use GK::Fold instead (as clearly separating folding info is tricky)

# last update. 190527 along with revision of msa_conv.pl


no warnings 'substr';

our $OKAY_FORMAT = { #1, both import and export (standard) // 2, export only , txt format, // 3 export only, fst format // 4, export one seq per file // 5, mask
	'fst'=>1,
	'nxs'=>1,
	'aln'=>1,
	'sto'=>1,
	'phy'=>1,
	'mase'=>1,
	'namelist'=>2, #txt
	'paml'=>2, #txt
	'pure'=>3, #fst
	'single'=>4, #fst, multi-file in one dir
	'mask'=>5, #has its own export subroutine
};

sub new {
	my ($class)=shift @_;
	my $self= bless {},$class;
	$self->{env}{_lastid}=0;
	$self->{env}{_maxlen}=0; # literal length with gaps / all characters
	$self->{env}{_maxlen_name}=0; # literal length with gaps / all characters
	return $self;
}
sub import { # import one file to Obj
	my ($self, $filepath,$plain) = @_; #NOTE. use $plain with caution
	if (!$filepath or !-e $filepath or -B $filepath) {
		return 0; #"can't properly open file or file doesn't exist";
	}
	open (my $fh, $filepath);
	my $format=guess_msa($filepath);
	if ($plain) {
		my @d0=File::Spec->splitdir($filepath);
		my $sname=pop @$0;
		$sname=~s/\.\S+$//g;
		$self->read_fh_plain($fh, $sname);
	}
	elsif ($format=~/fst/) { $self->read_fh_fst($fh); }
	elsif ($format =~/aln|sto/) { $self->read_fh_aln($fh); }
	elsif ($format =~/mase/) { $self->read_fh_mase($fh); }
	elsif ($format =~/phy/) { $self->read_fh_phy($fh); }
	elsif ($format =~/nxs/) { $self->read_fh_nxs($fh); }
	else { return 0; }
	$self->{env}{_format}=$format;
	$self->_maxlen_batch;
	$self->_maxnamelen_batch;
	$self->_clean_seq;
	return 1;
}
# below, output to files, no more format check! input file root name must be complete path
sub export { #standard formats: fst,aln,phy,nxs,sto,mase
	my ($self,$format,$param)=@_; #$usefullname for phy/sto/aln only
	my $ofilert=$param->{fileroot}||''; # file path: no extension if NOT fixfilename
	my $fixoutname=$param->{fixfilename}?1:0;
	my $usefullname=$param->{fullname}?1:0; #only for phy and aln format. sequence title character number
	my $ofilestd=$fixoutname?$ofilert:$ofilert.'.'.$format;
	my $fh2;
	open ($fh2, ">", $ofilestd);
	$fh2=*STDOUT if !$fh2;
	# --------- write all seqs ------------
	$format=lc $format;
	if ($format=~/fst/) { $self->write_fh_fst($fh2); }
	elsif ($format=~/aln/) { $self->write_fh_aln($fh2,$usefullname); }
	elsif ($format=~/sto/) { $self->write_fh_sto($fh2,$usefullname); }
	elsif ($format=~/nxs/) { $self->write_fh_nxs($fh2); }
	elsif ($format=~/mase/) { $self->write_fh_mase($fh2); }
	elsif ($format=~/phy/) {
		$self->write_fh_phy($fh2,$usefullname);
		# my $ofile3=$ofilestd.'_fullname.phy';
		# my $fh3;
		# open ($fh3, ">", $ofile3);
		# $fh3=*STDOUT if !$fh3;
		# $self->write_fh_phy($fh3,1); #a version with full name. as trim to first 10 chars isn't enough
		# $ofilestd .= ' | '. $ofile3;
	}
	else { return 'unknown format to convert'; }
	return $ofilestd;
}
sub export_stripgaps { #return output filename
	my ($self,$ofileroot,$all_nogap,$count_unclear,$fixoutname)=@_;
	my $ofile=$fixoutname?$ofileroot:(sprintf '%s.nogap_%s.fst', $ofileroot, ($all_nogap?'all':'least1'));
	my $fh2;
	open ($fh2, ">", $ofile);
	$fh2=*STDOUT if !$fh2;
	my $nogapcols=$self->getnogapcols(!$all_nogap,$count_unclear);
	my $wmask=$self->write_fh_mask($fh2,$nogapcols);
	if ($wmask) {
		return $ofile;
	} else {
		close ($fh2);
		unlink $ofile;
		return '!!not all columns are gap-free?';
	}
}
sub export_mask { #return output filename
	my ($self,$ofileprefix,$maskdata,$fixoutname)=@_; #maskdata=[name, (char1 , char2, ...])
	my $mname=shift @$maskdata;
	my $chars=_getmaskchars($maskdata);
	my $ofile=$fixoutname?$ofileprefix:(sprintf '%s.mask_%s%s.fst', $ofileprefix, $mname,($chars?'_'.(join "", sort @$chars):''));
	my $fh2;
	open ($fh2, ">", $ofile);
	$fh2=*STDOUT if !$fh2;
	my $plist=$self->getmaskcols($mname, $chars);
	my $wmask=$self->write_fh_mask($fh2,$plist);
	if ($wmask) {
		return $ofile;
	} else {
		close ($fh2);
		unlink $ofile;
		return '!!specified mask name and/or chars are not found';
	}
}
sub export_special { #other special formats, to single file: namelist, pure, paml
	my ($self,$fileroot,$format,$fixoutname)=@_;
	my $feval=$self->format_type($format);
	my $ofile=$fixoutname?$fileroot:(sprintf '%s.%s.%s',$fileroot,$format,($feval==2?'txt':'fst'));
	my $fh2;
	open ($fh2, ">", $ofile);
	$fh2=*STDOUT if !$fh2;
	if ($format=~/pure/) { $self->write_fh_pure($fh2); }
	elsif ($format=~/namelist/) { $self->write_fh_namelist($fh2); }
	elsif ($format=~/paml/) { $self->write_fh_paml($fh2); }
	else { return 'unknown output format'; }
	# phase=>\&write_fh_phase,
	return $ofile;
}
sub export_single { # one seq per file, remove all gaps
	my ($self,$rootfile,$fixoutname)=@_; #$fixoutname =1, if rootfile is already a dir
	my $odir;
	if ($fixoutname) {
		$odir=$rootfile;
	} else {
		my @rpath=File::Spec->splitpath($rootfile);
		my $fname=pop @rpath;
		$odir=File::Spec->catdir(@rpath, sprintf '%s_single', $fname);
	}
	if (!-d $odir) {
		mkpath $odir;
	}
	if (-d $odir) { #still no dir after trying to make it
		foreach my $id (1..$self->lastid) {
			next if !$self->contains($id);
			my $ofilename=sprintf '%0*d__%s.fst', (length $self->lastid), $id, safename($self->getname($id));
			my $ofile = File::Spec->catfile($odir, $ofilename);
			open (my $fh3, ">", $ofile);
			$self->_write_fh_pure_1seq($fh3, $id);
			close ($fh3);
		}
		return $odir;
	} else {
		return "can't access output dir";
	}
}

# --------- access by id ------
sub drop { #delete seq by id
	my ($self,$id)=@_;
	if ($self->contains($id)) {
		delete $self->{$id};
	}
}
sub getseq {
	my ($self,$id)=@_;
	if ($self->contains($id)) {
		return $self->{$id}{seq};
	} else {
		return '';
	}
}
sub getvna { #vna isn't analysed in this module.
	my ($self,$id)=@_;
	if ($self->contains($id)) {
		return $self->{$id}{vna}||'';
	} else {
		return '';
	}
}
sub getname {
	my ($self,$id)=@_;
	if ($self->contains($id)) {
		return $self->{$id}{name};
	} else {
		return '';
	}
}
sub contains { #contains this id
	my ($self, $id)=@_;
	if (defined $self->{$id}{seq}) {
		return 1;
	} else {
		return 0;
	}
}
sub has_vna {
	my ($self,$id)=@_;
	if ($self->contains($id) and $self->{$id}{vna}) {
		return 1;
	} else {
		return 0;
	}
}

# ------- access ENV
sub format { #returns the format of the LAST imported file if multiple are imported
	my $self=shift;
	return $self->{env}{_format}||'unknown';
}
sub lastid {
	my ($self,$raise,$forceid)=@_;
	if ($forceid and $forceid=~/\d+/ and $forceid>$self->{env}{_lastid}) { #for safe
		$self->{env}{_lastid}=$forceid;
	} elsif ($raise) {
		$self->{env}{_lastid}++;
	}
	return $self->{env}{_lastid}||0;
}
sub maxlen {
	my ($self,$newlen)=@_; #update is for the seq-id to be checked
	if ($newlen and $newlen > $self->{env}{_maxlen}) {
		$self->{env}{_maxlen}=$newlen;
	}
	return $self->{env}{_maxlen}||0;
}
sub maxlen_name {
	my ($self,$newlen)=@_; #update is for the seq-id to be checked
	if ($newlen and $newlen > $self->{env}{_maxlen_name}) {
		$self->{env}{_maxlen_name}=$newlen;
	}
	return $self->{env}{_maxlen_name}||0;
}
sub totalseq { #for phy,nxs. return actual written seq total numbre
	my ($self)=@_;
	return ((scalar keys %{$self}) -1);
}

# ------ format eval -------
sub chk_format { #to unitify writing for fst,nxs
	my ($self,$ofmt)=@_;
	$ofmt = '' if !$ofmt;
	$ofmt = lc $ofmt;
	$ofmt=~s/\s|\W//g;
	if ($ofmt=~/^(fst|fsa|fas|fasta)$/) { return 'fst'; }
	elsif ($ofmt=~/^(nxs|nexus)$/) { return 'nxs'; }
	# elsif ($ofmt=~/^(aln|sto|phy|mase)$/) { return $ofmt; } #formats
	# elsif ($ofmt=~/^(namelist|mask|pure|single|paml)$/) { return $ofmt; } #special output
	else { return $ofmt; }
}
sub guess_msa { #non OO // feed in filepath // target formats: fst, aln, phy, nxs, sto, mase
	my $filepath=shift;
	open (my $fh, $filepath) or return 0;
	my $fmt='';
	my $firstline=0;
	while (<$fh>) {
		next if /^\s*$/; #ignore empty lines, then raise line count. debatable
		$firstline++;
		if ($firstline==1) { # some formats use the first line as meta-header(?)
			if (/^\s*  \d+ \s+ \d+  \s*$/x) { $fmt='phy'; } #seq count, seq length
			elsif (/^\s*CLUSTAL/i) { $fmt= 'aln'; }
			elsif (/^#\s*STOCKHOLM/i) { $fmt= 'sto'; }
			elsif (/^#\s*nexus/i) { $fmt= 'nxs'; }
			elsif (/^;;/i) { $fmt= 'mase'; }
			elsif (/^>/) { $fmt='fst'; }
		}
		last if $fmt;
	}
	close ($fh);
	return $fmt;
}
sub format_type {
	my ($self,$format)=@_;
	$format=$self->chk_format($format)||'';
	return $OKAY_FORMAT->{$format}||0;
}

# ------ internal use ------------
sub _maxnamelen_batch {
	my ($self, $m,$n)=@_;
	$m=1 if !$m;
	$n=$self->lastid if !$n;
	foreach my $j ($m..$n) {
		next if !$self->contains($j);
		$self->maxlen_name(length $self->getname($j));
	}
}
sub _maxlen_batch { #get the max len of many seqs
	my ($self, $m,$n)=@_;
	$m=1 if !$m;
	$n=$self->lastid if !$n;
	foreach my $j ($m..$n) {
		next if !$self->contains($j);
		$self->maxlen(length $self->getseq($j));
	}
}
sub _namespace { #for shorten name string, for phylip format
	my ($self,$name,$writename,$namelen,$namespace)=@_;
	if ($namelen and $namelen eq '-1') {
		$namelen=$self->maxlen_name;
	}
	if (!$namelen or $namelen < 10) {
		$namelen=10;
	}
	$namespace=3 if !$namespace;
	my $namefield;
	if (!$writename) {
		$namefield=sprintf "%s", ' ' x ($namelen+$namespace);
	} else {
		if ((length $name) > $namelen) {
			$namefield = sprintf "%s%s", (substr $name, 0, $namelen), (' ' x $namespace);
		} else {
			$namefield = sprintf "%-*s", ($namelen + $namespace), $name;
		}
	}
	return $namefield;
}
sub safename {
	my ($name)=@_;
	$name=~s/\W/_/g;
	return $name;
}
sub _rmspace { #literal whitespaces
	my $seq=shift;
	$seq=~s/\s//g;
	return $seq;
}
sub _rmgaps { #of a sequence. gap=non-AtoZ,0to9,*
	my $seq=shift;
	$seq=~s/[^a-z0-9*]//ig;
	return $seq;
}
sub _clean_seq { # all ~to -
	my $m=shift;
	foreach my $id (1..$m->lastid) {
		next if !$m->contains($id);
		my $s=$m->getseq($id);
		$s=~s/~/-/g;
		$m->{$id}{seq}=$s;
	}
}

# ------- change nt format to RNA or DNA -------
# will do all lines, including annotation
sub to_dna {
	my ($m, $asr) = @_;
	my $c=GenKode::CentralDogma->new;
	foreach my $id (1..$m->lastid) {
		next if !$m->contains($id);
		my $s=$m->getseq($id);
		if (!$asr) { #NOT as rna
			$s=$c->r2dna($s);
		} else {
			$s=$c->d2rna($s);
		}
		$m->{$id}{seq}=$s;
	}
}
sub to_rna {
	my $m=shift;
	$m->to_dna(1);
}

# -------- block process ---------
sub rename { # case INSENSITIVE
	my ($self,$renamelist)=@_; #H ref. {'curr name'=>'new name'}.
	if (!$renamelist) {
		return 0;
	}

	#case insenstivie, so lc all keys in list
	foreach my $kk (keys %$renamelist) {
		$renamelist->{lc $kk}=$renamelist->{$kk};
	}
	my $renamed=0;
	foreach my $id (1..$self->lastid) {
		next if !$self->contains($id);
		if ($renamelist->{lc $self->getname($id)}) {
			$self->{$id}{name}=$renamelist->{$self->getname($id)};
			$renamed=1;
		}
	}
	return $renamed;
}
sub reorder { # case INSENSITIVE
	my ($self,$order)=@_; #A ref
	if (!$order) {
		return 0;
	}
	
	my $order2;
	my $max=scalar @$order;
	for (my $i=0; $i<$max; $i++) {
		my $j=$i+1; #id num as in $self, begin from 1
		# also lc all elem in $order
		my $el = lc $order->[$i];
		if (!$order2->{$el}) {
			$order2->{$el}=$j;
		}
	}
	my $lastid2=0;
	my $base=$self->lastid * 10; # x10 to create a gap between old set and new set. in case old ones aren't deleted properly
	my $reordered=0;
	foreach my $id (1..$self->lastid) {
		next if !$self->contains($id);
		my $tmp=dclone $self->{$id};
		$self->drop($id);
		my $newid=0;
		if ($order2->{lc $tmp->{name}}) { #assign new order
			$newid=$order2->{lc $tmp->{name}}+$base;
			$reordered=1;
		} else {
			$newid=$id+$base*100;
		}
		# printf "[%s] %s - [%i]\n", $id, $tmp->{name}, $newid;
		$self->{$newid}=dclone $tmp;
		$lastid2=$newid if $lastid2<$newid;
	}
	$self->lastid(0, $lastid2);
	return $reordered;
}
sub _getmaskchars { # non-OO
	my ($chars)=@_;
	my $chars2;
	foreach my $cc (@$chars) {
		push @$chars2, (split '', $cc);
	}
	return $chars2;
}
sub getmaskcols { #return A ref for positions shown in mask. mask name, case INSENSITIVE
	my ($self,$maskname, $chars)=@_;
	my $plist; # A ref, begin from 0 for easy substr
	my $mseq;
	#find the mask seq
	foreach my $id (1..$self->lastid) {
		next if !$self->contains($id);
		if (lc $self->getname($id) eq lc $maskname) {
			$mseq=$self->getseq($id);
			last;
		}
	}
	if ($mseq) {
		my $p=0;
		my $chars2=_getmaskchars($chars);
		while (my $cc=substr $mseq, $p, 1) {
			if ($chars2) {
				foreach my $char (@$chars2) {
					if (lc $char eq lc $cc) { #this char in mask line is needed
						push @$plist, $p;
						last;
					}
				}
			} else {
				if ($cc!~/[-\._\s~\,]/) { #not common chars for a gap
					push @$plist, $p;
				}
			}
			$p++;
		}
	}

	return $plist;
}
sub getnogapcols { #return A ref for all or atleast1 no-gap column positions, all seq in that col must be non-gap. if flagged, return col positions that at least one seq in that col has to be non-gap.
	my ($self,$atleast_1_nogap,$count_unclear,$is_nt)=@_;
	my $min=$atleast_1_nogap?1:$self->totalseq; #number of seq must have non-gap
	my $cols; #A ref
	for my $col (0.. ($self->maxlen-1)) {
		my $count=0;
		for my $id (1..$self->lastid) {
			next if !$self->contains($id);
			my $c=substr $self->getseq($id), $col, 1;
			if (!_isgap($c,$count_unclear,$is_nt)) { #this char isn't a gap
				$count++;
			}
			if ($count>=$min) {
				push @$cols, $col;
				last;
			}
		}
	}
	return $cols;
}
sub _isgap { #a guess.
	my ($c,$count_unclear,$is_nt)=@_; #for $count_unclear.. "N" and "X" for nt and aa are NOT treated as gap here.
	return 1 if !$c;
	if ($c=~/^[\s\._\-\W\d]$/) { #surely a gap. space, dot, underscore,dash,equal sign,non-word,digit
		return 1;
	} else {
		if ($count_unclear) {
			if ($is_nt and $c=~/n/i) { #count "N"
				return 1;
			} elsif (!$is_nt and $c=~/x/i) { #count "X"
				return 1;
			}
		}
	}
	return 0;
}
sub map_case { #input another seq (<= target, MUST NOT contain non-words ), convert case according to this "another seq".
	my ($self,$id,$refseq)=@_;
	my $s0=$self->getseq($id);
	return 0 if !$refseq;
	$refseq=~s/\W//g;
	$s0=~s/\W//g;
	if ($s0=~/$refseq/i) {
		my $pre=$`;
		my $s1='';
		my $last=(length $self->getseq($id))-1;
		my $cutpre=0;
		my $j=0;
		for my $i (0..$last) {
			my $c=substr $self->getseq($id), $i, 1;
			my $added=0;
			if ($c=~/\w/) {
				$cutpre++;
				if ($cutpre>(length $pre)) {
					if ($j<(length $refseq)) {
						my $cr=substr $refseq,$j,1;
						$added=1;
						if ($cr=~/[a-z]/) {
							$s1.=lc $c;
						} else {
							$s1.=uc $c;
						}
					}
					$j++;
				}
			}
			if (!$added) {
				$s1.=$c; #no change
			}
		}
		$self->{$id}{seq}=$s1;
	}
}


# -------- read files ------------
sub read_fh_plain { # all text in file is from one sequence chunk. ONLY remove white spaces, not taking responsibility to weird characters
	my ($self,$fh,$name) = @_;
	my $id=$self->lastid;
	$id++;
	my $seq='';
	while (<$fh>) {
		my $s=$_;
		$s=~s/\s//g;
		$seq.=$s;
	}
	$self->{$id}{name}=$name||'seq_'.$id;
	$self->{$id}{seq}=$seq;
	$self->lastid(0,$id);
	1;
}
sub read_fh_fst { #for standard fst and fasta-like nxs
	my ($self,$fh,$nxsflag) = @_;
	my $idbegin=$self->lastid;
	my $id=$idbegin;
	while (<$fh>) {
		chomp;
		next if /^\s*$/;
		if ($nxsflag and /^(END)?;/i) {
			last;
		}
		elsif (/^>/ or ($nxsflag and /^\[\d+\]\s+/)) {
			$id++;
			$self->{$id}{name}=$'||'seq_'.$id;
			$self->{$id}{seq}='';
		} else {
			$self->{$id}{seq} .= $_;
		}
	}
	# last record
	$self->lastid(0,$id);
	# $self->_maxlen_batch(($idbegin+1),$self->lastid);
	1;
}
sub read_fh_aln { #for aln, sto and aln-line nxs
# name seq...
# name2 seq2...
# 
# name seq_cont...
# name2 seq_cont2...
	my ($self,$fh,$nxs) = @_;
	my $idbegin = $self->lastid;
	my $id=$idbegin;
	while (<$fh>) {
		chomp;
		last if m{^//$}; #do not exist in clustal and nxs. safe
		last if /^;$/ and $nxs;
		next if /^\s*CLUSTAL|^#\s*STOCKHOLM/i;
		next if /^#/;
		if (/^\s*$/) {
			$self->lastid(0,$id); # make sure obj's lastid is right
			$id=$idbegin;
		}
		elsif (/^(\S+)\s+(\S+)\s*$/) {
			my ($name,$seq)=($1,$2);
			$id++;
			if (!$self->{$id}) {
				$self->{$id}{name}=$name;
				$self->{$id}{seq} = $seq;
			} else {
				$self->{$id}{seq} .= $seq;
			}
		}
	}
	# update maxlen
	# $self->_maxlen_batch(($idbegin+1),$self->lastid);
	1;
}
sub read_fh_phy {
# name seq...
# name2 seq2...
# 
#      seq_cont...
#      seq_cont2...
	my ($self,$fh) = @_;
	my $idbegin=$self->lastid;
	my $id=$idbegin;
	my $idlast=$idbegin;
	while (<$fh>) {
		chomp;
		next if /^\s*\d+\s+\d+$/; #1st line

		# reinitialise $id while reaching an empty line
		if (/^\s*$/) {
			$idlast=$id;
			$id=$idbegin;
		}
		else { $id++; }

		if (/^(\S+)\s+(\S.*\S)\s*$/) { # 1st block . with names
			my ($name,$seq)=($1,$2);
			$seq=_rmspace($seq);
			if (!$self->{$id}) {
				$self->{$id}{name}=$name;
				$self->{$id}{seq}=$seq;
			}
		}
		elsif (/^\s+(\S.*\S)\s*$/) { # following block, no names
			my $seq2=$1;
			$seq2=_rmspace($seq2);
			$self->{$id}{seq} .= $seq2;
		}
	}
	$self->lastid(0,$idlast);
	# $self->_maxlen_batch(($idbegin+1),$self->lastid);
	1;
}
sub read_fh_mase { # will add extra mask line if sites info exists
	my ($self,$fh) = @_;
	my $site; #A ref.
	my $siteflag; #for sites
	my $idbegin=$self->lastid;
	my $id=$idbegin;
	my $nameflag;
	while (<$fh>) {
		next if !/\S/;
		chomp;
		if (/^;no comment/) { #start of a new seq
			$id++;
			$siteflag=0;
			$nameflag=1; #name is always the next line
		} elsif ($siteflag and !/^;; \d+,\d+/) {
			$siteflag=0;
		} elsif (/^;;# of segments=\d+ (.+)$/ and !$siteflag) {
			$siteflag=1;
			push @$site,{ name=>$1, seq=>'' }; #sites
		} elsif (/^;;( \d+,\d+.+)$/ and $siteflag) {
			$site->[-1]{seq}.= $1;
		} elsif (/^;/) { #comments, do nothing
			next;
		} elsif ($nameflag) {
			$self->{$id}{name}=$_;
			$self->{$id}{seq}="";
			$nameflag=0;
		} else {
			$self->{$id}{seq}.=$_;
		}
	}
	#mk mask line if exists
	if ($site) {
		foreach my $sr (@$site) {
			$id++;
			my $map;
			my $i=1; #human order
			while ($sr->{seq}=~/(\d+),(\d+)/g) {
				while ($i<$1) { $map.='-'; $i++; }
				while ($i<=$2 and $i>=$1) { $map.='x'; $i++; }
				next if $i>$2;
			}
			$self->{$id}{name}=sprintf "MASK_%s",$sr->{name};
			$self->{$id}{seq}=$map;
		}
	}
	#unify
	$self->lastid(0,$id);
	# $self->_maxlen_batch($idbegin,$self->lastid);
	1;
}
sub read_fh_nxs {#contains ONLY seq blocks. Do NOT feed in complicated files with trees/etc
	my ($self,$fh) = @_;
	my $matrix;
	my $interleave;
	while (<$fh>) {
		chomp;
		next if /^\s*$/;
		next if /^#\s*nexus/i;
		if (/interleaved\s*=\s*yes/i or /interleave/i) {
			$interleave=1;
		} elsif (/^END;/i) {
			last;
		} elsif (/^matrix/i) { #start to intake data
			if ($interleave) { #process like ALN
				$self->read_fh_aln($fh,1);
			} else { #process like fst
				$self->read_fh_fst($fh,1);
			}
		}
	}
	1;
}

# ----------- write formats ----------
sub write_fh_fst {
	my ($self,$fh)=@_;
	foreach my $id (1..$self->lastid) {
		next if !$self->contains($id);
		$self->write_fh_fst_1seq($fh,$self->getname($id),$self->getseq($id));
		if ($self->getvna($id)) {
			$self->write_fh_fst_1seq($fh,$self->getname($id).'__vna',$self->getvna($id));
		}
	}
	1;
}
sub write_fh_fst_1seq {
	my ($self,$fh,$n,$s)=@_;
	printf $fh ">%s\n%s\n",$n,$s;
	1;
}
sub write_fh_phy {
#name max limit: 10 chars + 3 spaces
	my ($self,$fh, $fullname)=@_;

	# header line.
	printf $fh " %s %s\n",$self->totalseq,$self->maxlen;

	my $block=0;
	my $rem=$self->maxlen % 60;
	my $firstline=1; #phy only needs titles for the 1st block
	my $minus=60;
	while ($block < $self->maxlen) {
		if (($self->maxlen - $block)<=$rem) {
			$minus=$rem; #for end of file that's less than 60 chars
		}
		foreach my $id (1..$self->lastid) {
			next if !$self->contains($id);

			#process name for 1st block
			my $namefield=$self->_namespace($self->getname($id),$firstline,($fullname?0:-1), 5);
			# if ($fullname) {
				# my ($namelen,$namespace);
				# $namefield=$self->_namespace($self->getname($id),$firstline,$self->maxlen_name, 5);
			# } else {
				# $namefield=$self->_namespace($self->getname($id),$firstline);
			# }

			#process seq line
			my $thisline;
			$thisline=substr $self->getseq($id),$block,$minus;
			$thisline=~s/\s/-/g; #fill in "-" if no enough seq
			$thisline=~s/(.{10})/$1 /g; #10 chars per sub-block
			$thisline=~s/\s*$//g; #rm tailing space				

			#output
			printf $fh "%s%s\n",$namefield,$thisline;
		}
		printf $fh "\n";
		$firstline=0;
		$block+=60;
	}
	1;
}
sub write_fh_aln {
	my ($self,$fh,$usefullname)=@_;
	printf $fh "CLUSTALW multiple sequence alignment\n\n";
	$self->_write_alnblock($fh,$usefullname);
	1;
}
sub write_fh_sto {
	my ($self,$fh,$usefullname)=@_;
	printf $fh "# STOCKHOLM 1.0\n\n";
	$self->_write_alnblock($fh,$usefullname);
	printf $fh "//\n";
	1;
}
sub _write_alnblock { #for both aln/sto. WITHOUT header
#name max length =10. space=3
	my ($self,$fh,$usefullname,$namespace)=@_;

	my $block=0;
	my $rem=$self->maxlen % 60;
	my $minus=60;
	while ($block < $self->maxlen) {
		if (($self->maxlen - $block)<=$rem) {
			$minus=$rem; #for end of file that's less than 60 chars
		}
		foreach my $id (1..$self->lastid) {
			next if !$self->contains($id);
			#process name
			my $namefield=$self->_namespace($self->getname($id),1,($usefullname?-1:0),$namespace);
			#process seq line
			my $thisline='';
			$thisline=substr $self->getseq($id),$block,$minus;
			$thisline=~s/\s/-/g; #fill in "-" if no enough seq

			#output
			printf $fh "%s%s\n",$namefield,$thisline;
		}
		printf $fh "\n";
		$block+=60;
	}
	1;
}
sub write_fh_mase {
	my ($self,$fh)=@_;
	printf $fh ";; written by %s\n",__PACKAGE__;
	foreach my $id (1..$self->lastid) {
		next if !$self->contains($id);
		printf $fh ";no comment\n%s\n%s\n", $self->getname($id), $self->getseq($id);
	}
	1;
}
sub write_fh_nxs { #not interleaved // $param->{type} REQUIRED! << otherwise it may cause problem
# DataType = { standard | DNA | RNA | nucleotide | protein | continuous }
	my ($self,$fh,$type)=@_;
	my $typewarn=1;
	if (!$type or $type!~/standard | DNA | RNA | nucleotide | protein | continuous/xi) {
		$type='standard';
		$typewarn=2;
	}
	printf $fh "#NEXUS\n";
	my $time= localtime(time);
	printf $fh "[saved by %s on %s]\n",__PACKAGE__, $time;
	printf $fh "BEGIN DATA;\n  DIMENSIONS NTAX=%s NCHAR=%s;\n  FORMAT DATATYPE=%s\n  GAP=-\n  ;\n",
		$self->totalseq,
		$self->maxlen,
		uc $type;
	printf $fh "MATRIX\n";
	my $i=1;
	foreach my $id (1..$self->lastid) {
		next if !$self->contains($id);
		my $seq2;
		if ($self->getseq($id)=~/\[|\]/) { # if seq has "[" and "]", replace them!
			$seq2=$self->getseq($id);
			$self->getseq($id)=~tr/\[\]/()/;
		}
		printf $fh "[%i] %s\n%s\n",$i,$self->getname($id),$self->getseq($id);
		$i++;
	}
	printf $fh ";\nEND;\n";
	return $typewarn;
}

# -------- write special format ---------
sub write_fh_namelist { #<name><\t><id>
	my ($self,$fh)=@_;
	foreach my $id (1..$self->lastid) {
		next if !$self->contains($id);
		printf $fh "%s\t%03d\n",$self->getname($id)||'seq_'.$id,$id;
	}
	1;
}
sub write_fh_pure { #out put to fst only
	my ($self,$fh)=@_;
	foreach my $id (1..$self->lastid) {
		next if !$self->contains($id);
		$self->_write_fh_pure_1seq($fh, $id);
	}
}
sub _write_fh_pure_1seq {
	my ($self,$fh,$id)=@_;
	my $s2=_rmgaps($self->getseq($id));
	$self->write_fh_fst_1seq($fh,$self->getname($id),$s2);
}
sub write_fh_mask { #output to fst only. mask name is case insensitive
	my ($self,$fh,$plist)=@_;
	if (!$plist or scalar @{$plist} == 0) { #mask not found or no chars retained
		return 0;
	} else {
		foreach my $id (1..$self->lastid) {
			next if !$self->contains($id);
			my $tseq='';
			foreach my $p (@$plist) {
				if ($p<(length $self->getseq($id))) {
					$tseq.=substr $self->getseq($id), $p, 1;
				} else {
					$tseq.='.';
				}
			}
			$self->write_fh_fst_1seq($fh,$self->getname($id),$tseq);
		}
		return 1;
	}
}
sub write_fh_paml { #for paml's input sequence file. make sure input doesn't contain unnecessary seqs
	my ($self,$fh)=@_;
	printf $fh "%s %s\n",$self->totalseq,$self->maxlen;
	$self->write_fh_fst($fh);
}


# ------------ consensus -------------
sub mk_consensus {
	my ($self,$cutoff,$mode,$exclude)=@_; # cutoff (0,1], $exclude=H ref. case insensitive
	#mode: 0: basic, atcgn / 1: intermediate , atcgnYR /2: full. all symbols. only apply if >=$cutoff. if cutoff is too low, e.g. 25%, a case of A 70%  G 28% and T 2% will result [AG]=>[R]/[N].
	#but usually should work well if cutoff is decent (>=50%) since only one type of nt will be in this range
	foreach my $k (keys %$exclude) {
		my $k2=lc $k;
		if ($k2 ne $k) {
			$exclude->{$k2}=1;
			delete $exclude->{$k};
		}
	}
	if ($cutoff!~/^\d\.?\d*$/) {
		$cutoff=0.75; #default
	}

	my $consen='';
	foreach my $i (0.. ($self->maxlen-1)) {
		my $colstat=$self->get_column($i,$exclude);
		my $addnt=colstat_2_consen($colstat,$cutoff,$mode);
		$consen.=$addnt;
	}
	my $id2=$self->lastid+1;
	$self->{$id2}{name}=sprintf 'consensus_%i', ($cutoff*100);
	$self->{$id2}{seq}=$consen;
	$self->lastid(0,$id2);
}
sub get_column { # return H ref of all NT's from one col. {$id}->{nt}. MAKE SURE vna is NOT individual sequences // SKIPed seq's will be ignored // $pos is MACHINE ORDER
	my ($self,$pos,$exclude)=@_;
	my $col={
		gap=>0,
		a=>0,
		t=>0,
		c=>0,
		g=>0,
		u=>0,
		n=>0,
	};
	foreach my $s (1..$self->lastid) {
		next if !$self->contains($s);
		next if $exclude->{lc $self->getname($s)};
		my $char=(substr $self->getseq($s),$pos,1)||'-';
		$char=lc $char;
		if (!defined $col->{$char}) {
			$col->{gap}++;
		}
		elsif ($char!~/[atcgu]/i) { #including N and R,Y,B.....
			$col->{n}++;
		}
		else {
			$col->{$char}++;
		}
	}
	return $col;
}
sub colstat_2_consen {
	my ($colstat,$cutoff,$mode)=@_;
	my $a=$colstat->{a}||0;
	my $t=($colstat->{t}||0)+($colstat->{u}||0);
	my $c=$colstat->{c}||0;
	my $g=$colstat->{g}||0;
	my $n=$colstat->{n}||0;
	my $gap=$colstat->{gap}||0;
	my $total=$a+$t+$c+$g+$n+$gap;
	my $useu=$colstat->{u}?1:0;

	$cutoff=0.75 if !$cutoff;
	$mode=0 if !$mode;
	if ($cutoff<=0.5 and !$mode) {
		$mode=1; #lower cutoff will cause two or more nt have same ratio. in this case, simple mode won't work
	}

	my $consen='';
	if ($gap>($total-$gap)) {
		$consen='-';
	} else {
		#1. chk simple case
		my $set={
			a=>$a,
			t=>$t,
			c=>$c,
			g=>$g,
			n=>$n,
		};
		foreach my $this (keys %$set) {
			if (($set->{$this}/$total)>=$cutoff) {
				$consen=uc $this;
				last;
			}
		}
		#2.
		if (!$consen) {
			if (!$mode) {
				$consen='N';
			}
			else {
				if ($mode==1) {#YR
					my $yratio=($t+$c)/$total;
					my $rratio=($c+$a)/$total;
					if ($yratio>=$cutoff and $rratio<$cutoff) {
						$consen='Y';
					}
					elsif ($rratio>=$cutoff and $yratio<$cutoff) {
						$consen='R';
					}
					else {
						$consen='N';
					}
				}
				else { #all other symbols B,H,S....
					my $ratio;
					if ($a and $t) { $ratio->{w}=($a+$t)/$total; }
					if ($c and $g) { $ratio->{s}=($c+$g)/$total; }
					if ($a and $c) { $ratio->{m}=($a+$c)/$total; }
					if ($t and $g) { $ratio->{k}=($t+$g)/$total; }
					if ($t and $c and $g) { $ratio->{b}=($t+$c+$g)/$total; }
					if ($a and $t and $g) { $ratio->{d}=($a+$t+$g)/$total; }
					if ($a and $t and $c) { $ratio->{h}=($a+$t+$c)/$total; }
					if ($a and $c and $g) { $ratio->{v}=($a+$c+$g)/$total; }
					if ($a and $t and $c and $g) { $ratio->{n}=($a+$c+$g+$t+$n)/$total };
					foreach my $sym (keys %$ratio) {
						if ($ratio->{$sym}<$cutoff) {
							delete $ratio->{$sym};
						} else {
							$consen=uc $sym;
						}
					}
					if (!$consen or (scalar keys %$ratio)>1) {
						$consen='N';
					}
				}
			}
		}
	}
	if ($useu and $consen eq 'T') {
		$consen='U';
	}
	return $consen;
}


1;





__END__

#write phase
=pod
sub write_fh_phase { #will lose any "vna" info. only output a template (includes species # and length), need manual edit later
	my ($self,$fh,$param)=@_;
	#param->{type} : DNA, PROT CODON STRUCT
	#put skip/excluded in one section, that don't count in total seq number
	#put non-skipping seq in main section, counts in total seq number
	printf $fh "# raw template generated by %s\n",__PACKAGE__;
	printf $fh "# base sequence block for phase, remember to manually edit before use in program!\n\n";
	my ($main,$skip);
	foreach my $id (@{$self->get_ids}) {
		if ($self->{$id}{skip}) {
			push @$skip,$id;
		} else {
			push @$main,$id;
		}
	}
	if ($skip) {
		printf $fh "# skipped line block\n";
		printf $fh "# REMEMBER TO:\n# UNCOMMENT and MOVE sequences as required by PHASE!\n\n";
		foreach my $id (@{$skip}) {
			printf $fh "#%s\n#%s\n\n",$self->{$id}{name},$self->{$id}{seq};
		}
	}
	if ($main) {
		printf $fh "# main sequence block\n\n";
		printf $fh "%s %s %s\n\n",$self->_get_totalseq(1,1),$self->{0}{MAXLEN},uc $param->{type};
		foreach my $id (@{$main}) {
			printf $fh "%s\n%s\n\n",$self->{$id}{name},$self->{$id}{seq};
		}
	}
	1;
}

=cut
