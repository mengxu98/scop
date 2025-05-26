#!/bin/bash

# Exit on error
set -e

# Default value for re_run_pkgdown
RE_RUN_PKGDOWN=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --re-run-pkgdown)
            RE_RUN_PKGDOWN=true
            shift
            ;;
        *)
            echo "Unknown parameter: $1"
            echo "Usage: $0 [--re-run-pkgdown]"
            exit 1
            ;;
    esac
done

echo "Starting pkgdown site deployment to gh-pages..."

# Function to check if working tree is clean
check_git_clean() {
    if [[ $(git status --porcelain) ]]; then
        return 1
    else
        return 0
    fi
}

# Check if docs directory exists
if [[ ! -d "docs" ]]; then
    echo "docs/ directory not found. Will run pkgdown::build_site()..."
    RE_RUN_PKGDOWN=true
fi

# Save current branch name
CURRENT_BRANCH=$(git symbolic-ref --short HEAD)

# Check if there are uncommitted changes
if ! check_git_clean; then
    echo "Stashing uncommitted changes..."
    git stash push -m "Temporary stash before deploying pkgdown"
    CHANGES_STASHED=true
else
    CHANGES_STASHED=false
fi

# Function to cleanup and restore original state
cleanup() {
    echo "Cleaning up..."
    git checkout "$CURRENT_BRANCH"
    if [[ "$CHANGES_STASHED" == true ]]; then
        echo "Restoring stashed changes..."
        git stash pop
    fi
}

# Set trap to ensure cleanup runs on script exit
trap cleanup EXIT

# Build the site using pkgdown if requested or if docs doesn't exist
if [[ "$RE_RUN_PKGDOWN" == true ]]; then
    echo "Building pkgdown site..."
    Rscript -e 'pkgdown::build_site()'
else
    echo "Skipping pkgdown site build (use --re-run-pkgdown to force rebuild)"
fi

# Create and switch to a new temporary branch
echo "Creating temporary branch..."
git checkout --orphan temp_gh_pages

# Remove everything except docs/
git rm -rf .
git clean -fxd
mv docs/* .
rm -rf docs/

# Add all files
echo "Adding built files..."
git add .

# Commit changes
echo "Committing changes..."
git commit -m "Build pkgdown site"

# Force push to gh-pages branch
echo "Pushing to gh-pages branch..."
git push -f origin temp_gh_pages:gh-pages

echo "Deployment completed successfully!"
echo "Please ensure GitHub Pages is configured to deploy from the gh-pages branch."
echo "Your site will be available at: https://[username].github.io/[repository]/"
