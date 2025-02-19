# Copy what is done by the macports bootstrap.sh used on github
# runners to move homebrew out of the way:
# https://github.com/macports/macports-ports/blob/master/.github/workflows/bootstrap.sh

echo "===> START DEBUG:  python3:  $(which python3)"

# Move directories to /opt/*-off
echo "Moving directories..."
sudo mkdir /opt/local-off /opt/homebrew-off
test ! -d /usr/local || /usr/bin/sudo /usr/bin/find /usr/local -mindepth 1 -maxdepth 1 -type d -print -exec /bin/mv {} /opt/local-off/ \;
test ! -d /opt/homebrew || /usr/bin/sudo /usr/bin/find /opt/homebrew -mindepth 1 -maxdepth 1 -type d -print -exec /bin/mv {} /opt/homebrew-off/ \;

# Unlink files
echo "Removing files..."
test ! -d /usr/local || /usr/bin/sudo /usr/bin/find /usr/local -mindepth 1 -maxdepth 1 -type f -print -delete
test ! -d /opt/homebrew || /usr/bin/sudo /usr/bin/find /opt/homebrew -mindepth 1 -maxdepth 1 -type f -print -delete

# Rehash to forget about the deleted files
hash -r

echo "===> END DEBUG:  python3:  $(which python3)  python3-config:  $(which python3-config)"
