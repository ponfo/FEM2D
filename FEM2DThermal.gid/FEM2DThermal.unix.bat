#!/bin/bash

rm $2/projectData.dat

cat <<EOT>> projectData.dat
$1
$3
EOT

$3/precomp $2/$1
$3/main $2/$1
