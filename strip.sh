#!/bin/sh

sed 's/^.*chr4//' $1 | tr -d ' ' | cut -c1-9
