#!/bin/bash

# To generate the token
# goto gitlab, Click your icon: Edit profile->Access tokens->Add new token
# generate token with permissions: read_repository, write_repository

# Variables
GITLAB_USER="gabrielw26"
GITLAB_TOKEN="aaaa" # Ensure this token has the necessary permissions
WORKING_DIR=~/tmp/mirror
SOURCE_GITLAB_REPO_URL="https://gitlab.fizyka.pw.edu.pl/wtools/wslda.git"
SOURCE_GITLAB_WIKI_URL="https://gitlab.fizyka.pw.edu.pl/wtools/wslda.wiki.git"
TARGET_GITLAB_REPO_URL="https://$GITLAB_USER:$GITLAB_TOKEN@gitlab.com/coldatoms/wslda.git"
TARGET_GITLAB_WIKI_URL="https://$GITLAB_USER:$GITLAB_TOKEN@gitlab.com/coldatoms/wslda.wiki.git"

# Create working directory
rm -rf $WORKING_DIR
mkdir -p $WORKING_DIR
cd $WORKING_DIR || exit

# Clone the source GitLab repository
if git clone --mirror $SOURCE_GITLAB_REPO_URL; then
  cd wslda.git || exit

  # Remove existing target GitLab remote if it exists
  git remote remove target_gitlab 2>/dev/null

  # Add and push to target GitLab remote
  git remote add target_gitlab $TARGET_GITLAB_REPO_URL
  git push --mirror target_gitlab
else
  echo "Failed to clone source GitLab repository."
  exit 1
fi

# Clone the source GitLab wiki repository
cd $WORKING_DIR || exit
if git clone --mirror $SOURCE_GITLAB_WIKI_URL; then
  cd wslda.wiki.git || exit

  # Remove existing target GitLab remote if it exists
  git remote remove target_gitlab 2>/dev/null

  # Add and push to target GitLab remote
  git remote add target_gitlab $TARGET_GITLAB_WIKI_URL
  git push --mirror target_gitlab
else
  echo "Failed to clone source GitLab wiki repository."
  exit 1
fi

echo "Mirroring completed successfully."
