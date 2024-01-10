#!/usr/bin/env bash
zcat -c $1 | paste - - - - | cut -f 2 | tr -d '\n' | wc -c
