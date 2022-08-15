# Info
This folder contains libraries. In most cases these are clones from W-Tools repo.

# Update of repos
```bash
git subtree pull --prefix lib/wdata https://gitlab.fizyka.pw.edu.pl/wtools/wdata master --squash
git subtree pull --prefix lib/wderiv https://gitlab.fizyka.pw.edu.pl/wtools/wderiv master --squash
git subtree pull --prefix lib/winterp https://gitlab.fizyka.pw.edu.pl/wtools/winterp master --squash
git subtree pull --prefix lib/wbox https://gitlab.fizyka.pw.edu.pl/wtools/wbox master --squash
```

# Adding new repo
Change `pull` into `add`.  
If you get error:
```
Working tree has modifications.  Cannot add.
```
make a new clone of the repo and repeat the operation. 
