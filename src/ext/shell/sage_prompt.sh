# -*- shell-script -*-
# to be sourced

__sage_describe_branch()
{
    local git_dir=$(git rev-parse --git-dir 2>/dev/null) 
    if [ ! $git_dir ]; then
	return
    fi

    local color=
    local op=
    if [ -d "$git_dir/.dotest" ];  then
        if [ -f "$git_dir/.dotest/rebasing" ];  then
            op="rebase"
	    color=red
        elif [ -f "$git_dir/.dotest/applying" ]; then
            op="am"
	    color=red
        else
            op="am/rebase"
	    color=red
        fi
    elif [ -f "$git_dir/.dotest-merge/interactive" ] ;  then
        op="rebase-i"
	color=red
    elif [ -d "$git_dir/.dotest-merge" ] ;  then
        op="rebase-m"
	color=red
    elif [ -f "$git_dir/MERGE_HEAD" ];  then
        op="merge"
	color=red
    elif [ -f "$git_dir/index.lock" ];  then
        op="locked"
	color=red
    elif [ -f "$git_dir/BISECT_LOG" ]; then
	op="bisect"
	color=red
    else
	git diff --no-ext-diff --quiet --exit-code
	if [ $? -ne 0 ]; then
	    op="modified"
	    color=red
	fi
    fi

    branch="$(git symbolic-ref --short HEAD)"
    if [ $? -ne 0 ]; then
	branch="<detached-$(git log --format='%h' -1)>"
	if [ ! $color ]; then 
	    color=red
	fi
    else
	if [ ! $color ]; then 
	    color=green
	fi
    fi
    
    local status
    if [ $op ]; then
	status="$branch|$op"
    else
	status="$branch"
    fi
    
    local ansi_red='\e[0;31m'
    local ansi_green='\e[0;32m'
    local ansi_blue='\e[1;34m'
    local ansi_clear='\e[0m'
    if [ $color == 'red' ]; then
	echo "$ansi_red{$status}$ansi_clear"
    elif [ $color == 'green' ]; then
	echo "$ansi_green($status)$ansi_clear"
    else
	echo "$status"
    fi
}

__sage_ps1()
{
    local ansi_bblack='\e[1;30m'
    local ansi_green='\e[0;32m'
    local ansi_bblue='\e[1;34m'
    local ansi_blue='\e[0;34m'
    local ansi_clear='\e[0m'
    local time24h="\t"
    local pathshort="\w"
    echo "$ansi_bblue""Sage:$ansi_clear"\
	"$ansi_bblack$HOSTNAME$ansi_clear $time24h" \
	"\$(echo -e \$(__sage_describe_branch)) $ansi_blue$pathshort$ansi_clear\n\$ "
}

