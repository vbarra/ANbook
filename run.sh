git add *
git commit -m "changes"
git push
cd ..
jupyter-book build  --all ANbook/ 
cd ANbook
ghp-import -n -p -f _build/html
