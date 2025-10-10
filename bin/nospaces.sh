#!/bin/bash

function nospaces {
	for OLDNAME in *
	do
    	NEWNAME=$(echo "${OLDNAME}" | tr -s ' ' '_')
	    if [ "${OLDNAME}" != "${NEWNAME}" ]
	    then
	        mv "${OLDNAME}" "${NEWNAME}"
	    fi
	done
}

function surf_dirs {
	cd $1
	nospaces
	for D in *
	do
		if [ -d "${D}" ]
		then
			surf_dirs $(echo $(pwd)/${D})
			cd ..
		fi
	done
}

surf_dirs $1