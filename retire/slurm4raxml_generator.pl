#!/usr/bin/perl
use strict;
use warnings;

### to generate scripts for making raxmlHPC trees
### my first time trying of using ARGV and flags! -w-
###
### by kiyoko kisaragi @ http://www.pocchong.de
### created: 2013-7-29
### last modified: 2013-11-15

use File::Spec;
use File::Copy;
use File::Path;
use Getopt::Long;
# use KodeMagie::Method_HandyTools;
# use KodeMagie::Genes::Method_MSA;

###TOF#################

my %models=(
	CAT => 1,
	MIX => 1,
	MIXI => 1,
	GAMMA => 1,
	GAMMAI => 1,
	CAT_GAMMA => 1,
	CAT_GAMMAI => 1,
);
my %matrixes=(
	DAYHOFF => 1,
	DCMUT => 1,
	JTT => 1,
	MTREV => 1,
	LG => 1,
	WAG => 1,
	RTREV => 1,
	CPREV => 1,
	VT => 1,
	BLOSUM62 => 1,
	MTMAM => 1,
);


my (@files,$list,$outdir,$model,$invar,$type,$bootstrap,$empfreq,$dups,$matrix,@outgp);
my ($procs,$mem,$runday,$email);

GetOptions(
	'files|f=s{1,}'=>\@files,
	'list=s{1}'=>\$list,
	'outdir|o=s'=>\$outdir,
	'outgp|g=s{1,}'=>\@outgp,
	'type|t=s'=>\$type, #nt or aa
	'model|m=s'=>\$model, #gamma,cat, etc...
	'matrix|x=s'=>\$matrix,
	'empfreq|q'=>\$empfreq,
	'invar|i'=>\$invar,
	'day|runday|d=f'=>\$runday,
	'bootstrap|b=f'=>\$bootstrap,
	# 'procs|p=f'=>\$procs,
	# 'memory|r=f'=>\$mem,
	'email|e=s'=>\$email,
	'slice|pt=f'=>\$dups,
);
if ($list) {
	open (my $fh,$list);
	while (<$fh>) {
		chomp;
		push @files,$_ if -e $_;
	}
}

if (!@files or !$outdir) {
die <<USAGE;
----------------------------------
  *-f FILE1 FILE2 ...
  *-o OUT_DIR
  *-t TYPE #nt or aa
   -m MODEL #CAT/MIX/GAMMA [default]/CAT_GAMMA
   -x MATRIX #DAYHOFF/DCMUT/JTT/MTREV/WAG/RTREV/CPREV/VT/BLOSUM62/MTMAM/LG [default]
   -b BOOTSTRAP
   -i #active Invariable region computation
   -q #use empirical base frequencies
   -g OUTGROUP1 OUTGROUP2 ...

   -s SPLIT_PART #e.g. run same data for 10 times -> to get 1000 bootstrap in total...
   -d RUN_DAYs
#   -e EMAIL

* are required
----------------------------------

USAGE
}

mkpath $outdir unless -d $outdir;
my $pbsdir=File::Spec->catfile($outdir,'_slurm');
mkpath $pbsdir;

foreach my $file (@files) {
	next if (!-e $file or -z $file);
	my @fdir=File::Spec->splitpath($file);
	my $fname=$fdir[2]; $fname=~s/\.(\S+?)$//;

	printf "- %s\n", $file;

	# mk phy
	if ($file !~/\.phy/i) {
		print "input file has to be phy !\n";
		next;
	}
	
	my ($params,$pbsenv);
	{
		my @tmp=File::Spec->splitpath($outdir);
		$params->{phy}=$fname.'.phy';
		$params->{name}=$fname;
		$params->{model}=proc_model({model=>$model,invar=>$invar,type=>$type,matrix=>$matrix,empfreq=>$empfreq});
		$params->{bootstrap}=($bootstrap?$bootstrap:100);
		$params->{outgroup}=join ",",@outgp if @outgp;
		
		$pbsenv->{time}=$runday;
		# $pbsenv->{procs}=$procs;
		# $pbsenv->{mem}=$mem;
		$pbsenv->{email}=$email;
		$pbsenv->{projdir}=$tmp[2];
		$pbsenv->{fnamedir}=$fname;
		$pbsenv->{name}=$params->{name};
	}
	
	#copy phy into each working dir
	my ($phy,$outfile,$wkdir,@rxcmd,@pbshead,@out);
	
	if ($dups) {
		for (1..$dups) {
			my $loop=sprintf "%02d", $_;
			$wkdir=File::Spec->catfile($outdir,$fname.'_pt'.$loop);
			mkpath $wkdir;
			$phy=File::Spec->catfile($wkdir,$params->{phy});
			copy $file, $phy;
			$pbsenv->{name}.='_pt'.$loop;
			$pbsenv->{fnamedir}=$params->{name};
			@rxcmd=mk_raxml($params, ($_==1?0:1));
			@pbshead=proc_pbs_env($pbsenv);
			@out=(@pbshead,'',@rxcmd);
			$outfile=File::Spec->catfile($pbsdir,$params->{name}.'.slurm');
			write_file_unix($outfile,\@out);
		}
	}
	else {
		$wkdir=File::Spec->catfile($outdir,$fname);
		mkpath $wkdir;
		$phy=File::Spec->catfile($wkdir,$params->{phy});
		copy $file,$phy;
		@rxcmd=mk_raxml($params);
		@pbshead=proc_pbs_env($pbsenv);
		@out=(@pbshead,'',@rxcmd);
		$outfile=File::Spec->catfile($pbsdir,$params->{name}.'.slurm');
		write_file_unix($outfile,\@out);
	}	
}


sub write_file_unix {
	my ($path,$chunk)=@_; #A ref
	local ($\,$/);
	open (my $fh,">",$path);
	foreach (@$chunk) {
	binmode $fh;
		print $fh $_,"\n";
	}
}
sub mk_raxml {
	my ($params,$bootstraponly)=@_;#H ref. keys=phy,name,model,bootstrap,outgroup
	my $rseed=int(rand(10000))+1; #for -p seed
	my $bseed=int(rand(10000))+1; #for bootstrap seed
	$params->{bootstrap}=100 if (!$params->{bootstrap} or $params->{bootstrap}<=0);

	my $tree0=$params->{name}.'_tree0';
	my $tree_bs=$params->{name}.'_bs'.$params->{bootstrap};
	my $tree_final=$params->{name}.'_tree_final';

	my $tree0_cmd=sprintf "raxmlHPC-HYBRID-AVX2 -s %s -n %s -m %s -f d -i 10 -p %d -T 2 -# 2", $params->{phy},$tree0,$params->{model},$rseed;
	my $tree_bs_cmd=sprintf "raxmlHPC-HYBRID-AVX2 -s %s -n %s -m %s -b %d -# %d -i 10 -p %d -T 2", $params->{phy},$tree_bs,$params->{model},$bseed,$params->{bootstrap},$rseed;
	my $tree_final_cmd=sprintf "raxmlHPC-HYBRID-AVX2 -s %s -n %s -m %s -t %s -z %s -f b -i 10 -p %d",$params->{phy},$tree_final,$params->{model}, ( sprintf 'RAxML_result.%s.RUN.0', $tree0 ), 'RAxML_bootstrap.'.$tree_bs, $rseed;

	if ($params->{outgroup}) {
		foreach ($tree0_cmd,$tree_bs_cmd,$tree_final_cmd) {
			$_.=sprintf " -o %s",$params->{outgroup};
		}
	}
	my @cmd;
	if ($bootstraponly) {
		push @cmd,$tree_bs_cmd;
	} else {
		push @cmd,$tree0_cmd,$tree_bs_cmd,$tree_final_cmd;
	}
	@cmd;
}
sub proc_pbs_env {
	 my ($params)=@_; #H ref, str. keys=time,procs,mem,email,projdir*,fnamedir*
	 my @pbsout;
	 push @pbsout,
		'#!/bin/bash',
		'#SBATCH --account=def-zimmerly',
		(sprintf '#SBATCH --time=%i-00:00:00', $params->{time}||7), #use defined day or 7 days by default
		'#SBATCH --mem=1GB',
		'#SBATCH --nodes=1',
		'#SBATCH --ntasks=8',
		'#SBATCH --cpus-per-task=2',
		# (sprintf "#SBATCH --mail-user= %s",($params->{email}?$params->{email}:'kosmoflips@gmail.com')),
		'#SBATCH --mail-type=ALL',
		(sprintf '#SBATCH --job-name=%s', $params->{name}||''),
		'#SBATCH --output=%x-%j.out',
		'',
		'cd ~',
		(sprintf "cd %s",$params->{projdir}),
		(sprintf "cd %s",$params->{fnamedir}),
		'',
		'module load nixpkgs/16.09  gcc/5.4.0  openmpi/2.0.2',
		'module load raxml/8.2.11',
		'';
	@pbsout;
}
sub proc_model {
	my $params=shift; #H ref --> keys=model,invar,type,matrix,empfreq
	my $model;
	if (!$params->{model} or !$models{uc$params->{model}}) { $model='GAMMA'; }
	else { $model=uc$params->{model}; }
	$model=$model.'I' if $params->{invar};
	
	if ($params->{type} eq 'nt') { return 'GTR'.$model; }
	elsif ($params->{type} eq 'aa') {
		my $matrix;
		if (!$params->{matrix} or !$matrixes{uc$params->{matrix}}) { $matrix='LG'; }
		else { $matrix=uc$params->{matrix}; }
		$matrix=$matrix.'F' if $params->{empfreq};
		return ('PROT'.$model.$matrix);
	}
	else { 0; }
}
