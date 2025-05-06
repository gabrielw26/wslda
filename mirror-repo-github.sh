#!/bin/bash

# To generate the token
# goto github, Click your icon: Settings->Developer Settings->Personal access tokens->Classic
# generate token with permissions: repo, write:packages

# Variables
GITHUB_USER="gabrielw26"
GITHUB_TOKEN="aaaa" # Ensure this token has the necessary permissions
WORKING_DIR=~/tmp/mirror
GITLAB_REPO_URL="https://gitlab.fizyka.pw.edu.pl/wtools/wslda.git"
GITLAB_WIKI_URL="https://gitlab.fizyka.pw.edu.pl/wtools/wslda.wiki.git"
GITHUB_REPO_URL="https://$GITHUB_USER:$GITHUB_TOKEN@github.com/gabrielw26/wslda.git"
GITHUB_WIKI_URL="https://$GITHUB_USER:$GITHUB_TOKEN@github.com/gabrielw26/wslda.wiki.git"

# Create working directory
rm -rf $WORKING_DIR
mkdir -p $WORKING_DIR
cd $WORKING_DIR || exit

# Clone the GitLab repository
if git clone --mirror $GITLAB_REPO_URL; then
  cd wslda.git || exit

  # Remove existing GitHub remote if it exists
  git remote remove github 2>/dev/null

  # Add and push to GitHub remote
  git remote add github $GITHUB_REPO_URL
  git push --mirror github
else
  echo "Failed to clone GitLab repository."
  exit 1
fi

# Clone the GitLab wiki repository
cd $WORKING_DIR || exit
if git clone --mirror $GITLAB_WIKI_URL; then
  cd wslda.wiki.git || exit

  # Remove existing GitHub remote if it exists
  git remote remove github 2>/dev/null

  # Add and push to GitHub remote
  git remote add github $GITHUB_WIKI_URL
  git push --mirror github
else
  echo "Failed to clone GitLab wiki repository."
  exit 1
fi

echo "Mirroring completed successfully."
