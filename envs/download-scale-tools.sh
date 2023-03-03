#!/bin/sh
srcDir=$(dirname -- "$0")
dir=${1:-"$(readlink -f "${srcDir}/../bin")"}
mkdir -p $dir
echo "Downloading to $dir"
curl http://scale.pub.s3.amazonaws.com/scaleTools/cicd/bcParser/master/230220-gf9e5e90/bc_parser -o $dir/bc_parser && chmod 755 $dir/bc_parser

