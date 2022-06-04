package GenKode::Tree;
use strict;
use warnings;
use Data::Dumper;
# PhyloTree stuff
# may have lots of errors since tree format is worsely unified than MSA
# by kiyo @ http://www.pocchong.de
# created: 2014-1-23
# last modified: 14-7-5

use base 'Exporter';
our @EXPORT = qw/
read_tree_file
file_2newick
file_2nexus
file_rename
file_rminfo
/;

our @EXPORT_OK = qw/

/;

our %EXPORT_TAGS = (
);


#feed in FILE (NEWICK except the convert one), to outfile's fh
sub read_tree_file { #input 1 NEWICK file, return A ref of 1-line tree strings
	my $file=shift;
	local $/=undef;
	open (my $fh,$file);
	my $chunk=<$fh>;
	$chunk=~s/\s//g;
	my @ts=split /;/,$chunk;
	@ts=map { $_ .=';' } @ts; #add ending ";" back
	return \@ts;
}
sub file_2nexus {
	my ($file,$fh)=@_; #all trees should have same name set
	my $trees=read_tree_file($file);
	my $index=index_names($trees->[0]);
	write_nexus_header($index,$fh);
	my $count=1;
	foreach my $tree (@$trees) {
	#rename to numbers
		printf $fh "   tree converted_%d = ",$count;
		tree_rename($tree,$fh,$index);
		print $fh "\n";
		$count++;
	}
	print $fh "end;\n";
}
sub file_2newick { #only this has input file being nexus
	my ($file,$fh2,$keepnum)=@_;#if $keepnum, use nexus numbers as names
	open (my $fh, $file);
	my $flag;
	my $index;
	while (<$fh>) {
		if (!$keepnum) {
			if (!$flag and /^\s*translate\s*$/i) {
				$flag=1;
			} elsif ($flag and /(\d+)\t([^,]+)\,?$/) {
				$index->{$1}=$2;
			} elsif ($flag and /^\s*;\s*$/) {
				$flag=0;
			}
		}
		if (/tree  \s+ \S+? \s* = \s* (   \(.+\);   ) \s* $/x) {
			my $tree=$1;
			if ($tree=~/\[&prob|\[&/) { #bayesian format, kinda unique so no worry about others
				$tree=~s!(\d)\[.+?\]!$1!g;
				my @piece=split /\[|\]/,$tree;
				$tree='';
				for my $i (0..$#piece) {
					if ($piece[$i]=~/&prob=(.+?),.+?/) {
						my $prob=$1;
						if ($prob=~/\de.?\d/i) { #is in sci.notation
							$piece[$i]=sprintf "%.3g", $prob;
						} else {
							$piece[$i]=$prob;
						}
					}
					$tree.=$piece[$i];
				}
			}
			tree_rename($tree,$fh2,$index);
			print $fh2 "\n";
		}
	}
	1;
}
sub file_rename {
	my ($file,$fh,$names)=@_;
	my $trees=read_tree_file($file);
	foreach my $tree (@$trees) {
		tree_rename($tree,$fh,$names);
		print $fh "\n";
	}
	1;
}
sub file_rminfo {
	my ($file,$fh,$opt)=@_;
	my $trees=read_tree_file($file);
	foreach my $tree (@$trees) {
		tree_rm_info($tree,$fh,$opt);
		print $fh "\n";
	}
	1;
}

#string manipulation
sub tree_rm_info { #to fh.
#opt keys: br|bt|name , name only control remove name
	my ($tree,$fh,$opt)=@_; #string, H ref
	# die Dumper $opt;
	my $brs=split_branch($tree);
	foreach my $leaf (@$brs) {
		if ($opt->{br}) {
			$leaf=leaf_rm_branch($leaf);
		}
		if ($opt->{bt}) {
			$leaf=leaf_rm_bootstrap($leaf, ($opt->{bt_min}||0));
		}
		if ($opt->{name}) {
			$leaf=leaf_rename($leaf,undef,1);
		}
		print $fh $leaf;
		print $fh "," unless $leaf=~/;$/;
	}
	1;
}
sub tree_rename { #to fh.
	my ($tree,$fh,$names)=@_; #$names = H ref { oldname => newname}
	my $brs=split_branch($tree);
	foreach my $leaf (@$brs) {
		$leaf=leaf_rename($leaf,$names); #if name isn't given, do nothing
		print $fh $leaf;
		print $fh "," unless $leaf=~/;$/;
	}
	1;
}
sub split_branch { #feed in 1 tree string, return A ref of branches
	my @pcs=split /,/,shift;
	$_=~s/^\s*?|\s*?$//sg foreach (@pcs);
	\@pcs;
}
sub leaf_rm_branch { #feed in 1 branch, remove branch info
	my $br=shift; #scalar
	$br=~s/:  \-? \d+  \.?  \d*  (e[\+\-]?\d+)?//igx;
	$br;
}
sub leaf_rm_bootstrap { #feed in 1 branch, remove bootstrap info.
	my $br=shift; #scalar
	my $min=shift;
	if (!$min) { #only remove values smaller than this
		$min=1000;
	}
	# my $val;
	# if ($br=~/\[  (  \-? \d+  \.?  \d*  (e[\+\-]?\d+)?  )  \]/igx) {
		# $val=$1;
	# }
	# elsif ($br=~/\)  (  \-? \d+  \.?  \d*  (e[\+\-]?\d+)?  )  /igx) {
		# $val=$1;
	# }

	# "[\d+]" format
	# $br=~s/\[  \-? \d+  \.?  \d*  (e[\+\-]?\d+)?    \]//igx;
	# ")\d+" format
	# $br=~s/\)  \-? \d+  \.?  \d*  (e[\+\-]?\d+)?  /)/igx;
	
	while ($br=~/\[  ( \-? \d+  \.?  \d*  (e[\+\-]?\d+)?  )  \]/igx) {
		my $val=$1;
		if ($val<$min) {
			$br=~s/\[  $val \]//ix;
		}
	}
	while ($br=~/\)  (  \-? \d+  \.?  \d*  (e[\+\-]?\d+)?  )  /igx) {
		my $val=$1;
		if ($val<$min) {
			$br=~s/\)  $val  /)/ix
		}
	}
	return $br;
}
sub leaf_rename { #feed in 1 branch, replace to new name (rm name if empty).
	my ($br,$names,$rmname)=@_;#scalar, H ref, if no new name, rm it
	my $name=leaf_getname($br);
	if ($rmname) {
		$br=~s/$name//;
	} elsif ($names and $names->{$name}) {
		my $n2=$names->{$name};
		$br=~s/$name/$n2/;
	}
	$br;
}
sub leaf_getname {
	my $leaf=shift;
	$leaf=~s/^\(+//g;
	$leaf=~s/:.+$//g;
	$leaf=~s/\).+$//g;
	$leaf=~s/\[.+$//g;
	return ($leaf?$leaf:undef);
	# my ($name)=$leaf=~/^  \(*  (  [\w\-]+  ) /x;
	# return ($name?$name:undef);
}
sub index_names { #mk list of node names from plain tree.  name => order
	my $brs=split_branch(shift);
	my $index;
	my $count=1;
	foreach (@$brs) {
		my $name=leaf_getname($_);
		$index->{$name}=$count;
		$count++;
	}
	$index;
}
sub write_nexus_header { #to fh, stop right before tree NAME = (,,,,,); ...
	my ($index,$fh)=@_;
	print $fh "#NEXUS\n\n";
	print $fh "begin taxa;\n";
	printf $fh "\tdimensions ntax=%d;\n",scalar keys %$index;
	print $fh "\ttaxlabels\n";
	foreach my $taxa (sort {$index->{$a} <=> $index->{$b}} keys %$index) {
		printf $fh "\t\t%s\n",$taxa;
	}
	print $fh "\t\t;\nend;\n";
	print $fh "begin trees;\n";
	print $fh "\ttranslate\n";
	foreach my $taxa (sort {$index->{$a} <=> $index->{$b}} keys %$index) {
		printf $fh "\t\t%d\t%s,\n",$index->{$taxa},$taxa;
	}
	1;
}

1;