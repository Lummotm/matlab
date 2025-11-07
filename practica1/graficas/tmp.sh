#!/usr/bin/env bash

files=$(find . -type f -name "*.png")

convert files ~/tmp/test.pdf
zathura ~/tmp/test.pdf
