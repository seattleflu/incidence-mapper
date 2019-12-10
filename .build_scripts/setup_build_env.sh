#!/usr/bin/env bash
docker login idm-docker-staging.packages.idmod.org/sfim_build_env:latest --username "${bamboo_UserArtifactory}" --password "${bamboo_PasswordArtifactory@Q}"
