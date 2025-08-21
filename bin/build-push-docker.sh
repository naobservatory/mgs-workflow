#!/bin/bash
# Script to build all images in the docker directory and push them to DockerHub
# Before use, make sure to authenticate with DockerHub using 'docker login' with the securebio username
# Usage: ./bin/build-push-docker.sh (must be used from mgs-workflow root directory)

set -euo pipefail

# Configuration
DOCKERHUB_ORG="securebio"

get_image_name() {
    local dockerfile="$1"
    local filename=$(basename "${dockerfile}")
    # Remove .Dockerfile extension and convert to lowercase
    local image_name="${filename%.Dockerfile}"
    echo "${image_name}" | tr '[:upper:]' '[:lower:]'
}

build_and_push() {
    local dockerfile="$1"
    local image_name=$(get_image_name "${dockerfile}")
    local full_tag="${DOCKERHUB_ORG}/${image_name}:latest"
    
    echo "Building ${image_name} from ${dockerfile}..."
    # cross-platform build (just in case one is accidentally running e.g. on Apple Silicon)
    # Currently amd64 only because as of July 2025 nao-rpkg base image is not available for arm64
    docker buildx build --platform linux/amd64 -f "${dockerfile}" -t "${image_name}:latest" .
    docker tag "${image_name}:latest" "${full_tag}"
    
    echo "Pushing ${full_tag}..."
    docker push "${full_tag}"
    
    echo "Successfully built and pushed ${full_tag}"
}

# Authenticate with DockerHub (assumes docker login has been run)
echo "Note: Make sure you're logged into DockerHub with 'docker login' with securebio username."

# Find and process all Dockerfiles in docker directory
echo "Discovering Dockerfiles in docker/ directory..."
for dockerfile in docker/*.Dockerfile; do
    if [ -f "${dockerfile}" ]; then
        echo "Processing ${dockerfile}"
        build_and_push "${dockerfile}"
        echo ""
    fi
done

echo "All images built and pushed successfully!"