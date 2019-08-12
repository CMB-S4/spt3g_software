#!/bin/sh

# Usage: getvers.sh <tree to get version info from> <output file>

set -e

exec 1>$2

cd $1

echo '# AUTO-GENERATED FILE: DO NOT EDIT'
echo

echo '# Contains version control information (revision IDs, paths, etc.)'
echo '# as of the last time make was run.'
echo

git_tree_modified()
{
	# This function stolen from FreeBSD's newvers.sh and is copyright
	# 2018 Ed Maste, under the same license as this repository.

        # git diff-index lists both files that are known to have changes as
        # well as those with metadata that does not match what is recorded in
        # git's internal state.  The latter case is indicated by an all-zero
        # destination file hash.

        local fifo

        fifo=$(mktemp -u)
        mkfifo -m 600 $fifo || exit 1
        git diff-index HEAD > $fifo &
        while read smode dmode ssha dsha status file; do
                if ! expr $dsha : '^00*$' >/dev/null; then
                        rm $fifo
                        return 0
                fi
                if ! git diff --quiet -- "${file}"; then
                        rm $fifo
                        return 0
                fi
        done < $fifo
        # No files with content differences.
        rm $fifo
        return 1
}

if [ -d .svn ]; then
	if (svn info --show-item url 1>/dev/null 2>/dev/null); then
		url=$(svn info --show-item url)
		relurl=$(svn info --show-item relative-url)
		baserev=$(svn info --show-item revision)
	else
		url=$(svn info | grep '^URL: ')
		url=${url#URL: }
		repobase=$(svn info | grep '^Repository Root: ')
		repobase=${repobase#Repository Root:}
		relurl=^$(echo $url | cut -c $(echo $repobase | wc -c )-)
		baserev=$(svn info | grep '^Revision: ')
		baserev=${baserev#Revision: }
	fi
	echo upstream_url=\"$url\"

	# Following redundant for SVN, but nice to have compatibility with
	# the git case below. Since it is just compat with git, "trunk" is
	# called "master" here.
	if echo $url | grep -q '/branches/'; then
		echo upstream_branch=\"$(echo $relurl | sed 's/.*\/branches\///g')\"
	elif echo $url | grep -q '/trunk'; then
		echo upstream_branch=\"master\"
	else
		echo upstream_branch=\"\"
	fi

	echo revision=\"$(svnversion)\"
	if (svnversion | grep -q M); then
		echo localdiffs=True
	else
		echo localdiffs=False
	fi
	if (echo $url | grep -q github.com); then
		echo gitrevision=\"$(svn pg --revprop -r $baserev git-commit)\"
	fi
	if echo $url | grep -q '/tags/'; then
		echo versionname=\"$(echo $relurl | sed 's/.*\/tags\///g')\"
	elif echo $url | grep -q '/releases/'; then
		echo versionname=\"$(echo $relurl | sed 's/.*\/releases\///g')\"
	else
		echo versionname=\"\"
	fi
elif [ -d .git ]; then
	if (git rev-parse --abbrev-ref --symbolic-full-name @{u}>/dev/null 2>/dev/null); then
		remote_branch=$(git rev-parse --abbrev-ref --symbolic-full-name @{u})
		if (git remote get-url 1>&2 2>/dev/null); then
			echo upstream_url=\"$(git remote get-url "$(echo $remote_branch | cut -d / -f 1)")\"
		else
			echo upstream_url=\"$(git config remote."$(echo $remote_branch | cut -d / -f 1)".url)\"
		fi
		echo upstream_branch=\"$(echo $remote_branch | cut -d / -f 2)\"
	else
		echo upstream_url=\"UNKNOWN VCS\"
		echo upstream_branch=\"UNKNOWN VCS\"
	fi
	echo revision=\"$(git rev-parse --verify HEAD 2>/dev/null)\"
	echo gitrevision=\"$(git rev-parse --verify HEAD 2>/dev/null)\"
	if git_tree_modified; then
		echo localdiffs=True
	else
		echo localdiffs=False
	fi
	echo versionname=\"$(git tag -l --points-at HEAD 2>/dev/null)\"
else
	echo upstream_url=\"UNKNOWN VCS\"
	echo upstream_branch=\"UNKNOWN VCS\"
	echo revision=\"UNKNOWN VCS\"
	echo gitrevision=\"UNKNOWN\"
	echo localdiffs=True
	echo versionname=\"UNKNOWN\"
fi

