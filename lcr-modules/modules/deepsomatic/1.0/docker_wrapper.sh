#!/bin/bash

# DeepSomatic Docker wrapper script
# This script handles Docker permissions and setup for the DeepSomatic module

set -euo pipefail

# Check if Docker is available
if ! command -v docker &> /dev/null; then
    echo "Error: Docker is not installed or not available in PATH" >&2
    exit 1
fi

# Check if Docker daemon is running
if ! docker info &> /dev/null; then
    echo "Error: Docker daemon is not running or accessible" >&2
    exit 1
fi

# Pull the DeepSomatic image if not already available
DOCKER_IMAGE="${DOCKER_IMAGE:-google/deepsomatic:1.9.0}"

if [[ "$(docker images -q ${DOCKER_IMAGE} 2> /dev/null)" == "" ]]; then
    echo "Pulling DeepSomatic Docker image: ${DOCKER_IMAGE}"
    docker pull "${DOCKER_IMAGE}"
fi

# Run the DeepSomatic command
exec docker "$@"