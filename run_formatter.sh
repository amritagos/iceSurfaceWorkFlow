#!/bin/sh
find . -regex '.*py' -exec black {} \;