package mkAnnoDB;

use strict;
use warnings;

use DBI;

sub connectdb {
	my $database=shift; #assume the path is okay (no extra folders to be created)
	my $driver   = "SQLite";
	my $dsn = "DBI:$driver:dbname=$database";
	my $dbh = DBI->connect($dsn, { RaiseError => 1 })
		or die $DBI::errstr;
	return ($dbh);
}


# --------------------------- gtf to sqlite ---------------------------

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