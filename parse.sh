#!/bin/sh

sed '/^#/d' $1 | sed 's/^.*QC+//' | tr -d ' '
