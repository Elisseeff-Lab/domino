# Designate desired directories and files
SRC_FOLDER_PATHS=(data-raw man R tests vignettes)
SRC_FILE_PATHS=(DESCRIPTION LICENSE NAMESPACE NEWS.md README.md inst/CITATION)

# Get files from directories
FILES=$(find "${SRC_FOLDER_PATHS[@]}" -type f)

# Add indicated individual files
for F in "${SRC_FILE_PATHS[@]}"
do
FILES+=" ${F}"
done

echo "${FILES[@]}"

# Pass these to github:
git config --global user.name 'GitHub Action'
git config --global user.email 'action@github.com'
# Fetch branches
git fetch
# Checkout target branch                         
git checkout $TARGET_BRANCH
# Copy files from master branch
git checkout master -- $FILES
# Commit to the repository (ignore if no changes)
git add -A
git diff-index --quiet HEAD ||  git commit -am "update files"
# Push to remote branch
git push origin $TARGET_BRANCH 