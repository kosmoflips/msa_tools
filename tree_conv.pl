#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

use lib ".";
use GenKode::Tree;

#no dir mode as tree files have diff. ext
my (@treefiles,@list,$help);
my $rminfo;
my $conv;
my $rename;
GetOptions(
	'files=s{1,}'=>\@treefiles,
	'lists=s{1,}'=>\@list,
	'remove'=>\$rminfo->{true},
		'branch'=>\$rminfo->{br},
		'bootstrap'=>\$rminfo->{bs},
		'btmin=s{1}' => \$rminfo->{bt_min},
		'name'=>\$rminfo->{name},
	'rename=s'=>\$rename->{map},
		'merge'=>\$rename->{merge},
	'convert'=>\$conv->{true},
		'simple|newick'=>\$conv->{nwk},
		'xus'=>\$conv->{nex},
		'notranslate'=>\$conv->{notrans},
	'help'=>\$help,
);

if (@list) {
	foreach (@list) {
		parse_file_list($_,\@treefiles);
	}
}

my $error='';
if ($help or !@treefiles) {
	$error="need input file";
}
elsif ($rminfo->{true}) {
	if ($rminfo->{br} or $rminfo->{bs} or $rminfo->{name}) {
	} else {
	use Data::Dumper;
	print Dumper $rminfo;exit;
		$error="need to specify what to remove";
	}
}
elsif ($rename->{map} and !(-e $rename->{map})) {
	$error="need rename map";
}
elsif ($conv->{true}) {
	if (($conv->{nwk} and $conv->{nex}) or (!$conv->{nwk} and !$conv->{nex})) {
		$error="need proper convert format";
	}
}
if ($error) { print STDERR<<FLAG;
----------------------------------
[-f TREEFILE1 TREEFILE2 ...]
#if a csv is give, names should be the same

[-convert] #use either [-x] or [-s]
  [-x] #convert to the complicated nexus format
  [-s] #convert to the simpler newick format.
    [-notr] #use with [-s]. use numbers from nexus as names directly

[-rename CSV_PATH] #CASE SENSITIVE. TAB separated. NO header. <CURR_NAME><\\t><NEW_NAME>
  [-m] #merge mode: new name will be APPENDED instead of replaced

[-remove]
  [-name] #remove all names
  [-branch] #remove branch lengths
  [-bootstrap] #remove bootstrap values
    [-btmin] #only remove values less than given number
      #e.g.: "75" means values from 0~74 will be removed

----------------------------------

- $error

FLAG
exit;
}

foreach my $t (@treefiles) {
	next if (!-e $t or -z $t);
	print ">$t. . .\n";
	if ($conv->{true}) {
		if ($conv->{nwk}) { #from neuxs to newick
			my $outfile=$t.'_newick.txt';
			open (my $fh2,">",$outfile);
			file_2newick($t,$fh2,$conv->{notrans});
		} else { #from newick to neuxs
			my $outfile=$t.'_nexus.txt';
			open (my $fh2,">",$outfile);
			file_2nexus($t,$fh2);
		}
	} elsif ($rename->{map}) {
		my $index=parse_name_idx($rename->{map},$rename->{merge});
		my $outfile=$t.'_renamed.txt';
		open (my $fh2,">",$outfile);
		file_rename($t,$fh2,$index);
	} elsif ($rminfo->{true}) {
		my $outfile=$t.'_rminfo.txt';
		open (my $fh2,">",$outfile);
		file_rminfo($t,$fh2,{br=>$rminfo->{br},bt=>$rminfo->{bs},name=>$rminfo->{name},bt_min=>($rminfo->{bt_min}||0)});
	}
}

sub parse_name_idx { #input tab separated file, output H ref. e.g. {oldname}=>{new name}
#you're responsible for the usibility of csv file!
	# shift if $_[0] eq __PACKAGE__;
	my ($csv,$merge)=@_;#file path, merge mode. <old><new>
	open (my $fh, $csv);
	my $names;
	while (<$fh>) {
		chomp;
		next if !/\S/;
		my @tmp=split /\t/;
		# if (!$merge) {
			# $names->{safe_name($tmp[0])}=safe_name($tmp[1]);
		# } else {
			# $names->{safe_name($tmp[0])}=safe_name($tmp[0]).'__'.safe_name($tmp[1]);
		# }
		if (!$merge) {
			$names->{$tmp[0]}=$tmp[1];
		} else {
			$names->{$tmp[0]}=($tmp[0] || '').'__'.($tmp[1] || '');
		}
	}
	$names;
}
