#!/bin/bash
docker build -t registry.gitlab.com/transipedia/dekupl-run:base -f Dockerfile.base .
docker push registry.gitlab.com/transipedia/dekupl-run:base
