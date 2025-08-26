#!/bin/bash
curdir=$(pwd)
for dir in autogen_md_input_files  autogen_resp_input  gm2md  gm  scoretraj  validation; do
	echo "Compiling ${dir}"
	path="../Internal/${dir}"
	cd ${path}
	./compile.sh
	cd ${curdir}
done
echo "Complete"
