package GTFsupport;

use strict;
use warnings;
use Data::Dumper;
use Storable qw/:DEFAULT nstore dclone/;
use DBI;

use base 'Exporter';
our @EXPORT = qw/
connectdb
get_exon_data_from_sqlite
get_exon_data_from_exons
map_exon_position
conv_coord_trx2gen
/;

our @EXPORT_OK = qw/

/;

our %EXPORT_TAGS = (
# translation=>[qw//]
);

# ----------- gtf extraction -------------
sub parse_gtf_attr {
	# need gene/trx id and exon number
	my ($attrline, $needed_attrs)=@_; # $needed_attrs, A ref of names to be returned
	#gene_id "ENSG00000284662"; gene_version "1"; transcript_id "ENST00000332831"; transcript_version "4"; exon_number "1"; gene_name "OR4F16"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "OR4F16-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS41221"; exon_id "ENSE00002324228"; exon_version "3"; tag "basic"; transcript_support_level "NA (assigned to previous version 3)"
	my $attr;
	if (!$needed_attrs) {
		$needed_attrs = [qw/gene_id  gene_version  transcript_id  transcript_version  exon_number/ ];
	}
	for my $item (@$needed_attrs) {
		if ($attrline=~/$item "(.+?)"/) {
			$attr->{$item}=$1;
		}
	}
	return $attr;
}

# ----------- coordinate conversion -------------
sub map_exon_position { # input genomic start/end for all exons, add transcriptomic start/end for each
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
sub conv_coordinate { # convert one position
	my ($exons, $pos, $is_trx_pos)=@_;
	foreach my $i (1..((scalar @$exons)-1)) {
		my $e=$exons->[$i];
		my ($p,$q); # original coord
		my ($m,$n); # map to coord
		# here assumes p,q,m,n are given in order, from small to big
		if ($is_trx_pos) { # trx to genomic
			($p,$q)=($e->[2],$e->[3]);
			($m,$n)=($e->[0],$e->[1]);
		} else { # genomic to trx
			($p,$q)=($e->[0],$e->[1]);
			($m,$n)=($e->[2],$e->[3]);
		}
		if ($pos>=$p and $pos<=$q) { # in this block
			if ($exons->[0]==1) { # plus strand
				return [$i, $m+($pos-$p)];
			} else {
				($p,$q)=($q,$p); # adjust for minus strand
				return [$i, $n-($pos-$q) ];
			}
		}
	}
	return [0, 0]; # position couldn't be converted, e.g. out of range due to wrong id, etc.
}
sub conv_coord_trx2gen {
	my ($p, $q, $exons)=@_; # p,q are transcript-based start/end, $exons has the format from &map_exon_position, may span over many exons
	if ($p>$q) {
		($p,$q)=($q,$p);
	}
	my $p2=conv_coordinate($exons, $p, 1);
	my $q2=conv_coordinate($exons, $q, 1);
	my $span;
	# die Dumper $exons;
	if ($p2->[0] < $q2->[0]) { # spans over many exons
		for my $j ($p2->[0]..$q2->[0]) {
			if ($p2->[0]==$j) { # starting exon, use converted p2 and exon end
				push @$span, [$j, $p2->[1], $exons->[$j][1]];
			}
			elsif ($q2->[0] == $j) {  # ending exon, use exon start and converted q2
				push @$span, [$j, $exons->[$j][0], $q2->[1]];
			}
			else { # middle exons
				push @$span, [$j, $exons->[$j][0], $exons->[$j][1]];
			}
		}
	} else { # on same exon
		push @$span, [$p2->[0], $p2->[1], $q2->[1]];
	}
	return $span;
}

# -------------- exon hash related -------------
sub get_exon_data_from_exons {
	my ($data, $trid, $gid)=@_;
	my $trxdata;
	if ($gid) {
		$trxdata=get_trx_data_by_gene($data, $gid, $trid);
	} else {
		$trxdata=get_trx_data_without_gene($data, $trid);
	}
	return $trxdata;
}
sub get_trx_data_by_gene { # find trx data by gid and trx id, no need to loop all gene ids
	my ($data, $gid, $tid)=@_;
	if ($data->{$gid}) {
		my $gdata=$data->{$gid};
		# find trx in current gene
		if ($gdata->{$tid}) {
			return get_trx_data($gdata->{$tid}, $data->{$gid}{info});
		}
	}
	return undef;
}
sub get_trx_data_without_gene { # find trx data by trx id only
	my ($data, $tid)=@_;
	foreach my $gid (keys %$data) {
		if ($data->{$gid}{$tid}) { # trx id exists in current gid
		# die Dumper $data->{$gid},999;
			return get_trx_data($data->{$gid}{$tid}, $data->{$gid}{info} );
		}
	}
	return undef;
}
sub get_trx_data { # used internally. get trx data from gene data
	my ($tdata, $info)=@_;
	$tdata->[0]={
		chr=>$info->{chr},
		strand=>$info->{strand} eq '1'?1:2,
		transcript_version=>$tdata->[0]{transcript_version}
	};
	return $tdata;
}


# -------------- db related -------------------
sub connectdb {
	my $database=shift; #assume the path is okay (no extra folders to be created)
	my $driver   = "SQLite";
	my $dsn = "DBI:$driver:dbname=$database";
	my $dbh = DBI->connect($dsn, { RaiseError => 1 })
		or die $DBI::errstr;
	return ($dbh);
}
sub get_exon_data_from_sqlite {
	my ($dbh, $trid, $DATA)=@_; # trid has no version
	# $DATA should be a global H-ref so any modification here will affect it globally
	if (!defined $DATA->{$trid}) {
		# update global $DATA
		get_exon_info_by_trx_id($trid, $DATA, $dbh);
		
		# trx id is found
		if ($DATA->{$trid}) {
			my $info=$DATA->{$trid};
			# die Dumper $info,23543;
			# SELECT gene_info.chr, exon.start, exon.end, exon.num, trx_info.id, trx_info.transcript_ver, gene_info.strand, gene_info.gene_id, trx_info.transcript_id, gene_info.gene_ver FROM ...
			my $info2;
			foreach my $eid (keys %$info) {
				if (!$info2) { # empty, add info into [0]
					push @$info2, {
						chr=>$info->{$eid}[0],
						transcript_version=>$info->{$eid}[5],
						strand=>$info->{$eid}[6]
					};
				}
				$info2->[$eid]=[$info->{$eid}[1],$info->{$eid}[2]];
			}
			# get exon info only
			$info2=map_exon_position($info2);
			$DATA->{$trid}=dclone $info2;
		}
	}
	if (!defined $DATA->{$trid}) { # this trx ID isn't found in db at all
		$DATA->{$trid}=undef;
		return ([0,[]]);
	} else {
		return $DATA->{$trid};
	}
}
sub get_exon_info_by_trx_id {
	my ($trx_name, $DATA, $dbh)=@_; # e.g. ENST0001 without version. $DATA is global
	my $sth=$dbh->prepare("SELECT gene_info.chr, exon.start, exon.end, exon.num, trx_info.id, trx_info.transcript_ver, gene_info.strand, gene_info.gene_id, trx_info.transcript_id, gene_info.gene_ver FROM exon JOIN trx_info ON exon.trx_id = trx_info.id JOIN gene_info ON gene_info.id = trx_info.gid WHERE trx_info.transcript_id = ?");
	$sth->execute($trx_name);
	while (my $tmp=$sth->fetch) {
		# die Dumper $tmp;
		$DATA->{$trx_name}{$tmp->[3]}=dclone $tmp; # key is exon-id
	}
	1;
}

# ---------------- gtf to sqlite ---------------------
sub gtf2db {
	my ($gtffile,$dbfile,$resume, $sampletest)=@_;
	if ($sampletest) { # sampletest : get 100 genes from GTF only, for testing purpose. overwrites resume option!
		$resume=0;
	}

	my $limit=100; # get sample only, 100 genes from gtf
	if ($sampletest) {
		$dbfile=$dbfile.'.samples.sqlite';
	} else {
		$limit=0;
	}

	if (-e $dbfile) {
		if ($limit) { #get samples, always re do from scratch
			unlink ($dbfile);
		}
		elsif (!$resume) { # remove existing file if not resuming unfinished building
			unlink ($dbfile);
		}
	}
	my $dbh=connectdb($dbfile);
	if (!$resume or $limit) { # make table only when db is new
		mktable($dbh);
	}
	# if errors occur it will die here.

	open (my $fh, $gtffile);
	my $currgid=0;
	my $currtid=0;
	my $resume_trx_id;
	if ($resume) {
		# find last transcript_id << because gtf file is parsed line by line, DO MAKE SURE it only works if input gtf is exactly the same as previous one.
		my $sth;
		$sth=$dbh->prepare("SELECT trx_info.id,transcript_id,gene_info.id from trx_info JOIN gene_info ON gene_info.id = trx_info.gid ORDER BY trx_info.id DESC LIMIT 1");
		$sth->execute();
		my ($last_trxid,$last_transcript_id, $last_gid)=@{$sth->fetch};

		$currgid=$last_gid; # so this gene won't be inserted again
		$resume_trx_id=$last_transcript_id;

		# remove all existing records for this id << easier than setting flags for missing ones
		$sth=$dbh->prepare("DELETE FROM exon WHERE trx_id=?");
		$sth->execute($last_trxid);
		$sth=$dbh->prepare("DELETE FROM cds WHERE trx_id=?");
		$sth->execute($last_trxid);
		$sth=$dbh->prepare("DELETE FROM trx_info WHERE id=?");
		$sth->execute($last_trxid);
	}
	my $resume_position_reached;
	if ($resume and $resume_trx_id) {
		print "    searching for resume point ...\n";
	}
	while (<$fh>) {
		next if /^#/;
		if ($resume and $resume_trx_id) { # need to resume AND know where to resume
			if (!$resume_position_reached) { # found the location in gtf to resume.
				if (/$resume_trx_id/) { # only do the string match before the location
					$resume_position_reached=1;
					printf "\n\n------- found resume point %s! -------\n\n", $resume_trx_id;
				} else {
					next;
				}
			}
		}
		my @c=split /\t/;
		my ($transcript_id,$transcript_ver)=$c[8]=~/transcript_id "(\S+)"; transcript_version "(\d+)"/;
		if ($c[2]=~/gene/) {
			my ($gene_id)=$c[8]=~/gene_id "(\S+)"/;
			printf "\n----%s----\n", $gene_id;
			if ($limit) { # sampling mode.
				$limit--;
				last if $limit==0;
			}
		}
		elsif ($c[2]=~/transcript/) {
			printf "    %s\n", $transcript_id;
			my ($gene_id,$gene_ver)=$c[8]=~/gene_id "(\S+)"; gene_version "(\d+)"/;
			my $sth;

			# see if this gene is saved and has $currgid
			$sth = $dbh->prepare("SELECT gene_id FROM gene_info WHERE id=?");
			$sth->execute($currgid);
			my $insertginfo=0;
			my $x0=$sth->fetch;
			if (!$x0 or scalar @$x0==0 or $x0->[0] ne $gene_id) {
				$insertginfo=1;
			} else {
				$insertginfo=0;
			}

			if ($insertginfo) {
				# insert gene info
				$sth = $dbh->prepare("INSERT INTO gene_info (gene_id,gene_ver,chr,start,end,strand) VALUES (?,?,?,?,?,?)");
				$sth->execute($gene_id, $gene_ver, $c[0], $c[3], $c[4], ($c[6]=~/\+/?1:2));

				# get gene table id
				$sth=$dbh->prepare("SELECT id FROM gene_info WHERE gene_id LIKE ?");
				$sth->execute($gene_id);
				$currgid=$sth->fetch->[0];
			}
			# insert trx info , can do this without checking trxID is because gtf file should be in order of transcript>exon/cds>etc...
			$sth=$dbh->prepare("INSERT INTO trx_info (gid, transcript_id, transcript_ver) VALUES (?,?,?)");
			$sth->execute($currgid, $transcript_id,$transcript_ver);

			# get trx table id
			$sth=$dbh->prepare("SELECT id FROM trx_info WHERE transcript_id LIKE ?");
			$sth->execute($transcript_id);
			$currtid=$sth->fetch->[0];
		}
		else {
			if ($c[2]=~/exon/) {
				my ($eid)=$_=~/exon_number "(\d+)"/;
				my $sth=$dbh->prepare("INSERT INTO exon (trx_id, num, start, end) VALUES (?,?,?,?)");
				$sth->execute($currtid,$eid,$c[3],$c[4]);
			}
			elsif ($c[2]=~/(start|stop)_codon|CDS/) {
				my $type; # 0 for cds, 1 for start, 2 for end
				if ($c[2]=~/CDS/) {
					$type=0;
				} elsif ($c[2]=~/start/) {
					$type=1;
				} else {
					$type=2;
				}
				my $frame; # use 0 for "type==end", use 1,2,3 for normal RF
				if ($type==0) { # CDS
					$frame=$c[7]+1;
				} else {
					$frame=0;
				}
				my $sth=$dbh->prepare("INSERT INTO cds (trx_id, type, frame, start, end) VALUES (?,?,?,?,?)");
				$sth->execute($currtid, $type, $frame, $c[3], $c[4]);
			}
		}
	}
	$dbh->disconnect;
	return 1;
}
sub mktable {
# table info
=pod
>>gene_info
- id*
- gene_id
- gene_ver
- chr
- start
- end
- strand

>>trx_info
- id*
- gid # link to gene_info.id
- transcript_id
- transcript_ver

>>exon
- id*
- trx_id # link to trx_info.id
- num
- start
- end

>>cds
- id*
- trx_id # link to trx_info.id
- frame # 1,2,3
- type
- start
- end
=cut
	my $dbh=shift;
	my $tb1_stat=qq(CREATE TABLE "gene_info" (
"id"	INTEGER NOT NULL UNIQUE,
"gene_id"	TEXT NOT NULL UNIQUE,
"gene_ver"	INTEGER,
"chr"	TEXT,
"start"	INTEGER NOT NULL,
"end"	INTEGER NOT NULL,
"strand"	INTEGER NOT NULL,
PRIMARY KEY("id" AUTOINCREMENT)
););
	my $tb1b_stat=qq(CREATE TABLE "trx_info" (
"id"	INTEGER NOT NULL UNIQUE,
"gid"	TEXT NOT NULL,
"transcript_id"	TEXT NOT NULL UNIQUE,
"transcript_ver"	INTEGER,
PRIMARY KEY("id" AUTOINCREMENT)
););
	my $tb2_stat=qq(CREATE TABLE "exon" (
"id"	INTEGER NOT NULL UNIQUE,
"trx_id"	TEXT NOT NULL,
"num"	INTEGER NOT NULL,
"start"	INTEGER NOT NULL,
"end"	INTEGER NOT NULL,
PRIMARY KEY("id" AUTOINCREMENT)
););
	my $tb3_stat=qq(CREATE TABLE "cds" (
"id"	INTEGER NOT NULL UNIQUE,
"trx_id"	TEXT NOT NULL,
"frame"	INTEGER NOT NULL,
"type"	INTEGER NOT NULL,
"start"	INTEGER NOT NULL,
"end"	INTEGER NOT NULL,
PRIMARY KEY("id" AUTOINCREMENT)
););

	my $rv1 = $dbh->do($tb1_stat);
	my $rv1b = $dbh->do($tb1b_stat);
	my $rv2 = $dbh->do($tb2_stat);
	my $rv3 = $dbh->do($tb3_stat);
	if($rv1<0 or $rv1b<0 or $rv2<0 or $rv3<0) {
	   die $DBI::errstr;
	} else {
	   return 1;
	}
}

1;
