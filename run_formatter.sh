#!/bin/sh
find ./scripts -regex '.*py' -exec black {} \;