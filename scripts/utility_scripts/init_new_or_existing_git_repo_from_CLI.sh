#!/usr/bin/bash

## sed -i 's/\r$//' init_new_or_existing_git_repo_from_CLI.sh

### to create a new git repo on the CLI

# change "new_project_repo" to my respective repo name
new_project_re="trial_init_command"


echo "# ${new_project_re}"
echo "# ${new_project_re}" >> README.md
# git init
# git add README.md
# git commit -m "first commit"
# git branch -M main
# git remote add origin https://github.com/HildaWacu/$new_project_repo.git
# git push -u origin main

### to push an exixsting repo from the CLI
# git remote add origin https://github.com/HildaWacu/$new_project_repo.git
# git branch -M main
# git push -u origin main