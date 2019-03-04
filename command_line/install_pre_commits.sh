#!/bin/bash
# LIBTBX_SET_DISPATCHER_NAME dev.dials.install_pre_commits

libtbx.pip install -U pip # can be removed with base moving to 2.7.16
libtbx.pip install pre-commit # can be removed when dependency is declared

for module_directory in $(libtbx.show_repository_paths); do
for module in $(libtbx.list_modules); do
if [ -d ${module_directory}/${module} ] && [ -e ${module_directory}/${module}/.pre-commit-config.yaml ]; then
  echo
  echo Installing pre-commits for repository ${module}
  pushd ${module_directory}/${module} >/dev/null
  libtbx.python -m pre_commit install
  sed -i '/import sys/aif not os.getenv("LIBTBX_BUILD"):\n  print("\\033[1;31mlibtbx environment not loaded. skipping pre-commits\\033[0m")\n  sys.exit(0)' .git/hooks/pre-commit
  popd >/dev/null
fi
done
done
