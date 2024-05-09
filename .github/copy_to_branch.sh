# Designate desired directories and files
FILES=(
'data-raw/'*
'man/'*
'R/'*
'tests/'*
'vignettes/'*
'DESCRIPTION'
'NAMESPACE'
'NEWS.md'
'README.md'
'inst/CITATION'
)

echo "${FILES[@]}"

# Pass these to github:
git config --global user.name 'GitHub Action'
git config --global user.email 'action@github.com'
# Fetch branches
git fetch
# Checkout target branch                         
git checkout $SRC_BRANCH
# copy files from the branch the action is being run upon
for F in ${FILES}; do
git checkout $SRC_BRANCH -- ${F}
done
# Commit to the repository (ignore if no changes)
git add -A
git diff-index --quiet HEAD ||  git commit -am "update files"
# Push to remote branch
git push origin $TARGET_BRANCH 