#!/bin/bash
echo "$DOCKER_PASSWORD" | docker login -u "$DOCKER_USERNAME" --password-stdin
docker build -t transipedia/dekupl-run .
docker push transipedia/dekupl-run