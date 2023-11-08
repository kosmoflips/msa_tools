# ////////////////////////////////////////////////////
# ----- commonly used non-graphic libs -----
# ////////////////////////////////////////////////////
library(stringr)
library(tidyr)
library(tidyverse)
library(dplyr)
library(reshape) # for "melt"
#library(gtools) # for combination


# ----- system/dir/file -----

# get current time
get_time <- function() { # get current time, formatted for mkdir use
	return(format(Sys.time(), "%y%m%d_%H%M%S"))
}

# save df to tsv
save_df_to_tsv <- function (data, file) {
	write.table(data, file = file, sep = "\t", quote = F, row.names = F )
}

# create path
mkpath <- function (path, mksubdir=0) {
	# if mksubdir is TRUE, will make a subdir under 'path' using current timestamp
	if (!dir.exists(path)) {
		dir.create(path)
	}
	
	if (mksubdir==1) {
		sdir=file.path(path, get_time())
		if (!dir.exists(sdir)) {
			dir.create(sdir)
		}
		return (sdir)
	}
	
	return(path)
}

# quick string match
string_contains <- function (string, pattern) { # tag, input string; pattern, for search
	if (is.na(str_match(string,pattern))==T ) { # no match
		return (FALSE)
	} else {
		return (TRUE)
	}
}
