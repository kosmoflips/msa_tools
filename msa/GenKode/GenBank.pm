package GenKode::GenBank;
use strict;
use warnings;

# for getting stuff from genbank

# http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
# id => UID list. Either a single UID or a comma-delimited list of UIDs may be provided. All of the UIDs must be from the database specified by db. There is no set maximum for the number of UIDs that can be passed to EFetch, but if more than about 200 UIDs are to be provided, the request should be made using the HTTP POST method.
# so I set up limit to be 100 here , for pubmed

#genbank record ref: http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html

use File::Spec;
use File::Path;
# use File::Temp;
use Data::Dumper;
use LWP::Simple;
# use XML::Simple;

#ref  ftp://ftp.ncbi.nlm.nih.gov/blast/documents/NEWXML/xml2.pdf

use Exporter 'import';

our @EXPORT = qw/
fetchseq
/;

our @EXPORT_OK = qw/

/;


sub fetchseq {#download a genbank entry into one given formats, RETURN downloaded file path or url if fail.
#full url: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=189430232&strand=1&seq_start=1&seq_stop=600&rettype=fasta&retmode=text');
	my ($param,$force)=@_; #acc (gi or acc),start,end,str,type,mode,db,odir,name

	if (!$param->{acc}) {
		print STDERR "need to specify id (gi or genbank accession)";
		return 0;
	}
	if ($param->{odir} and !-d $param->{odir}) {
		mkpath $param->{odir};
	}
	if (!$param->{odir} or !-d $param->{odir}) {
		print STDERR "need to specify output dir";
		return 0;
	}

	my $db=$param->{db}||'nuccore';
	my $type=$param->{type}||'fasta';
	my $mode=$param->{mode}||'text';
	my $p=$param->{start}||0;
	my $q=$param->{end}||0;
	($q,$p)=($p,$q) if ($p>$q);
	my $str=1;
	if ($param->{strand} and $param->{strand} ne '1') {
		$str=2;
	}

	my $url='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi';
	$url.="?db=".$db;
	$url.="&id=".$param->{acc};
	$url.="&strand=".$str;
	$url.="&seq_start=".$p if $p;
	$url.="&seq_stop=".$q if $q;
	$url.=sprintf "&rettype=%s&retmode=%s", $type,$mode;

	my $fname;
	if (!$param->{name}) {
		$fname=$param->{acc};
		$fname.='__p_'.$p if $p;
		$fname.='__q_'.$q if $q;
		$fname.='__str_'.$str if $str;
		$fname.='__'.$type if $type;
		$fname.='.'.$mode if $mode;
	} else {
		$fname=sprintf "%s_%s.%s", $param->{name},$type,($mode eq 'text'? 'txt':$param->{mode});
	}
	my $saveto=File::Spec->catfile($param->{odir}, $fname); #full path
	if (-e $saveto and !-z $saveto and !$force) { #force is used to overwrite existing file.
		return -1;
	}

	my $retry=0;
	while ($retry<=5) {
		if (is_success(getstore($url,$saveto))) {
			return $saveto; #success, end of sub
		} else {
			print STDERR " - failed... retry in 10 seconds...\n";
			sleep 10 unless $retry==5;
		}
		$retry++;
	}
	if ($retry) {
		print STDERR "have retried over 5 times, give up\n";
		unlink $saveto if -e $saveto;
		return 0;
	}
}

1;

__END__


###shared
sub xml2hash { #internal HERE
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
		return undef;
	} else {
		return $ref;
	}
}


####NCBI GenBank seq
# no more gi's as GB stops them in 2016-Sep or around
=pod
sub acc2gi {#convert ncbi accession number to gi (the very first gi appearing in the list)
#eg http://www.ncbi.nlm.nih.gov/nuccore/ABIY02000120.1?report=gilist&log$=seqview&format=text
	my ($acc,$db)=@_;
	$db='nuccore' if !$db; #works for both NT and AA while giving a confirmed gi number, I believe.
	my $url=sprintf "http://www.ncbi.nlm.nih.gov/%s/%s?report=gilist",$db,$acc;
	#above seems working fine. or NCBI's full url for this: http://www.ncbi.nlm.nih.gov/protein/EDU99216.1?report=gilist&log$=seqview&format=text
	my $saveto="$acc.txt";
	my $tmp=File::Temp->new(template=>'seq_XXXXX',suffix=>'.txt');
	close($tmp);
	my $giraw=getstore($url,$tmp->filename);
	local $/=undef;
	open (my $fh,$tmp->filename);
	my $chunk=<$fh>;
	close ($fh);
	my ($gi)=$chunk=~m|<pre>(\d+).*</pre>|s;
	return ($gi?$gi:undef);
}
=cut
# sub guess_gid { #a gi or an accession num? always return GI
	# my ($id)=@_;
	# if ($id=~/^\d+$/) { return $id; }
	# else { # ($id=~/^[a-z]+_?\d+(?:\.\d+)?$/i) 
		# return acc2gi($id);
	# }
# }

sub gb2href { #convert genbank xml to H ref
	xml2hash(shift,['GBSeq','GBSeqid','GBReference','GBAuthor','GBFeature','GBQualifier']);
}

### NCBI Blast Suite <== Local search
# require soft, e.g. ncbi-blast-2.2.29+-win64.exe being installed and accessable from cmd shell
# user manual: http://www.ncbi.nlm.nih.gov/books/NBK279690/
sub _verifyblastdb { #for local blast
	my ($db,$remote)=@_;
	if ($remote) {
		return -1 if !$db;
		my %BLAST_DB_REMOTE=( #COMMON DB only so far!
			nr=>1,
			nt=>1,
			refseq_rna=>1,
			refseq_protein=>1,
			refseq_genomic=>1,
			'bacteria and archaea'=>1,
			swissport=>1,
			chromosome=>1,
			pdb=>1,
			#no WGS as it needs more options
		);
		if (!$BLAST_DB_REMOTE{lc $db}) {
			return 0;
		}
	} else {
		return 0 if !$db;
		my @BLAST_DB_EXT=qw/.phr .pin .psq/;
		foreach (@BLAST_DB_EXT) {
			print "\n\n\n\n",$db.$_;die;
			if (!-e $db.$_) {
				warn "\n!! the given db misses one or more of *.phr,pin,psq files.\n $db";
				return 0;
			}
		}
	}
	return 1;
}
sub _verifyblastprog {
	my $method=shift;
	return 0 if !$method;
	my @BLAST_METHOD=qw/blastn blastp blastx tblastn/; #may add more in future
	foreach (@BLAST_METHOD) {
		return 1 if $method eq $_;
	}
	warn "\nthe given blast method isn't supported!\n";
	return 0;
}
sub localblast { #NOTE: should be ver. 2.2.31+ or higher, since I want xml2 format!
# if both subject (fasta file) and db are given. use subject first.
# other arg's not specified here can be given through ->{other}  as a string.

# http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
#input hash keys:: db||subject,prog,query,out
# when search remotely, use NR db only
#default output is XML for the ease of parse
#all arg's are same as ncbi's command line ver

# to specify a local db, use the path without "pal/etc" extension. make sure to match search program

	my $arg=shift; # db, query (fasta OR pure sequence), prog
	return undef if ref $arg ne "HASH";
	#required Hash keys

	return undef if (!$arg->{query} or !-e $arg->{query}); #must be a pure sequence file or fasta file!
	return undef if (!$arg->{out}); #path when save the blast result
	my $cmd=sprintf '-query "%s" -out "%s" -%s "%s" -outfmt %s%s',
		$arg->{query},
		$arg->{out},
		($arg->{subject}?'subject':'db'),
		($arg->{subject}?$arg->{subject}:$arg->{db}),
		($arg->{txtout}?0:14), #simple search. outformat=14 is xml2
		($arg->{remote}?' -remote':'');

	#optional
	if ( # task, for blastn only: 'blastn' 'blastn-short' 'dc-megablast' 'megablast' 'rmblastn'
		( ($arg->{prog} and $arg->{prog}=~/blastn/i) or
		($arg->{prog_path} and $arg->{prog_path}=~/blastn/i) ) and
		$arg->{task}
	) { $cmd .= sprintf " -task %s", $arg->{task}; }
	$cmd .= sprintf " -max_hsps %i", $arg->{max_hsps} if $arg->{max_hsps};
	$cmd .= sprintf " -max_target_seqs %i", $arg->{max_target_seqs} if $arg->{max_target_seqs};
	$cmd .= sprintf " -evalue %s", $arg->{evalue} if $arg->{evalue};
	$cmd.= sprintf ' %s ',$arg->{param} if $arg->{param};
# die $cmd;
	my @cmd;
	push @cmd,($arg->{prog_path} || $arg->{prog}),$cmd;
	eval { system(@cmd); };

	if ($@) { return $@; }
	else { 1; }
}
sub blast2href { #works for both xml=5 or 14 // read in XML blast result, return converted to hash ref. Feed in xml absolute path
	xml2hash(shift,['HitDescr','Hit','Hsp'], 1);
}

###Pfam
#pfam graphic: http://pfam.xfam.org/generate_graphic
#use CURL on Pfam: http://pfam.xfam.org/help#tabview=tab10
our $PSUB_STAT={ #return code from &pfam_submit, when not success
	1=>'curl error',
	2=>'server busy',
	3=>'invalid query',
	4=>'unknown error',
};
our $PRET_STAT={ #return code from &pfam_retrieve, when not success
	2 =>'xml + json',
	1 =>'xml',
	202 =>'server busy',
	0=>'fail',
};
sub pfam_submit {#read usage comments ↓ before use!!
#return: job id if success / 1, curl fail / 2, busy / 3, invalid seq / 4, other error
#submit local SEQUENCE-ONLY plain seq to pfam, return: job_id if success, err_code (self defined) if fail
#REQUIRES CURL BEING INSTALLED AND ACCESSABLE VIA CMD
	my ($seqfile,$outfile,$curl)=@_; #read_in, write_to , curl path on windows
	#my $searchurl='http://pfam.sanger.ac.uk/search/sequence';
	#link modified on 2014May2
	my $searchurl='http://pfam.xfam.org/search/sequence';
	my $arg=sprintf "%s -LH 'Expect:' -F seq='<%s' -F output=xml '%s' -o %s",($curl?$curl:'curl'),$seqfile,$searchurl,$outfile;
	if ($^O=~/mswin32/i) { #under HarveOS, use double quote instead of single
		$arg=~s/'/"/g;
	}
	my $i=0;
	my $returncode;
	while ($i<5) {
		system ($arg);
		#now analyse outcome
		if ($?) { #fail from CURL
			$returncode=1;
		} elsif (-z $outfile) { #empty outfile , something is wrong but now sure why
			$returncode=4;
		} else { #see if the submission (saved file) is good
			open (my $fh, $outfile);
			my $err;
			while (<$fh>) {
				chomp;
				next if /^\s*$/;
				$err=1 if /^<!DOCTYPE/;
				if ($err) {
					if (m{too many jobs}) { #server busy, retry later
						$returncode=2;
					} elsif (m{There was a problem determining the sequence type}) { #invalid letter in seq. ATCGN for dna, NO rna
						$returncode=3;
					}
				} else {
					if (m{<job job_id="(.+?)">}) {
						$returncode=$1;
					}
				}
				last if $returncode;
			}
		}
		$returncode= 4 if !$returncode;
		if ($returncode ne '4' and $returncode ne '2') { #no need to retry
			last;
		} else {
			$i++;
			last if $i<5;
			print "retrying. . . wait for 30 seconds. . .\n";
			sleep 30;
		}
	}
	return $returncode;
}
sub pfam_retrieve { #retrieve a job id to xml/json. protein only//requires JSON
#return: 2: json and xml; 1: xml only; 202: busy; 0: fail
	my ($id,$saveto,$savetojson)=@_;
	my $i=0;
	my $retry;
	my $url=sprintf "http://pfam.sanger.ac.uk/search/sequence/resultset/%s",$id;
	my $json=sprintf "http://pfam.sanger.ac.uk/search/sequence/graphic/%s",$id;
	my $stat;
	while ($i<5) {
		$stat=getstore($url, $saveto);
		if ($stat==200 and ((-s $saveto)>10)) { #retrieve success
			$stat=getstore($json,$savetojson);
			if ($stat>=200 and $stat<300) {
				return 2;
			} else {
				return 1;
			}
		} else {
			$i++;
			last if $i<5;
			print "retry. . . wait for 30 seconds. . .\n";
			sleep 30;
		}
	}
	#if code goes here, no good result, delete files if they exist
	unlink $saveto if -e $saveto;
	unlink $savetojson if -e $savetojson;
	return ($stat==202?$stat:0);
}
sub pfam_parse_json { #json to H ref
	my $file=shift; #pfam json file
	local $/ = undef;
	open (my $fh, $file);
	my $jref=decode_json(<$fh>); #from Cpan->JSON
	$jref;
=pod
	
	my ($region,$motif)=$chunk=~/
		"regions"\s*:\s*\[(.*?)\].*
		"motifs"\s*:\s*\[(.*?)\]
		/x; #dont need markup for now
		#print $region,$motif;die;
	my $info; #H ref
	$region.=$motif;
	my $count=1;
	while ($region=~/"metadata"\s*:\s*\{(.+?)\}/g) {
		my $tmp=$1;
#		print $tmp,"\n";die;
		$info->{$count}{text}=$1 if $tmp=~/"description"\s*:\s*"(.*?)"/;
		$info->{$count}{acc}=$1 if $tmp=~/"accession"\s*:\s*"(.*?)"/;
		$info->{$count}{id}=$1 if $tmp=~/"identifier"\s*:\s*"(.*?)"/;
		$info->{$count}{type}=$1 if $tmp=~/"type"\s*:\s*"(.+?)"/;
		$info->{$count}{start}=$1 if $tmp=~/"start"\s*:\s*"(.*?)"/;
		$info->{$count}{end}=$1 if $tmp=~/"end"\s*:\s*"(.*?)"/;
		$info->{$count}{db}=$1 if $tmp=~/"database"\s*:\s*"(.*?)"/;
		$info->{$count}{score}=$1 if $tmp=~/"score"\s*:\s*"(.*?)"/;
		#unify hash
		foreach (qw/text acc id type start end db score/) { $info->{$count}{$_}='-' if !$info->{$count}{$_};}
		$info->{$count}{text}=~s/\,/_/g;
		#method_peek($info);<>;
		$count++;
	}
	$info;
=cut
}
1;