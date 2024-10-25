#!/bin/sh
find . -type f -name "*.py" ! -name "*.npy" -exec black {} \;