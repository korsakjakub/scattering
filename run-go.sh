#!/bin/env zsh
docker run --rm -it -v $(pwd)/config.json:/go/src/config.json -v $(pwd)/go:/go/src/app scattering-go
