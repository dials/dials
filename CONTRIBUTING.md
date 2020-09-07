# Contributing to DIALS

We're happy to consider contributions from outside sources; whether in the form
of code, tickets, documentation or even just typo corrections! DIALS addresses
a wide range of use cases - so, if in the event that you are planning any large
scale work that you will want to merge back in, please contact us beforehand so
that we can discuss any potential impacts this might have.

Listed in this document are code standards and conventions that we try
to adhere to, some of which are essential and others that are just encouraged.
The intention is that all of the code should try to converge towards these.

## Configuring DIALS for development

1.  Install DIALS and its cctbx_project dependencies for development; For
    Linux/macOS, the current best way to create a fresh installation of DIALS
    and all of its dependencies is with the following commands:
    ```
    git clone https://github.com/dials/dials modules/dials
    python modules/dials/installer/bootstrap.py
    ```
2.  Activate the environment with `source <root>/dials`. This will
    need to be done every time you work on DIALS code.

The DIALS repository is now checked out in `<root>/modules/dials`. During
development, run tests with pytest with `libtbx.pytest --regression` to ensure
that all of the tests still pass.

If you update source code or change dependencies, you may occasionally need
to regenerate the static libtbx ecosystem and rebuild any C++ code. You can
go this by running `make reconf` in the `<root>/build` directory.

## Code Development Guidelines

These should be followed wherever possible - but remember the [Zen of Python]'s
"Practicality beats purity.", and [PEP8]'s "A Foolish Consistency is the
Hobgoblin of Little Minds" - it's okay to stray from the rules if you have a
good reason for it, but "I prefer it this way" isn't a strong enough reason -
there is real value to a standard style, because diverging from it sends a
message that the code in question is special and care should be taken.

- Err on the side of [PEP8] when making any style decision. In particular,
  **use PEP8 as a guide for naming** when you aren't sure the correct form to
  use - lowercase variable names, CamelCase class names etc.
- **Common imports should go at the top of the file**.  This makes it very easy
  to reason about the dependencies of a module, makes it hard to make
  conflicting definitions, and avoids duplicate imports scattered throughout
  the file. If an import is for an optional dependency, consider using a
  try/except block, with some fallback to identify the missing case (such as
  setting the name to None). Matplotlib is an exception to this guideline -
  because its startup logic defines the backend it uses, this can be imported
  inline. There are also a few exceptions to help avoid circular imports that
  are hard to remove. Further, unique imports are allowed within code branches
  or within optional functions/classes that are not always used when the file is
  imported.  This reduces the runtime import load for functionality that isn't
  universally used.
- Try not to do `from <module> import *` imports - it makes it hard to trace
  where definitions are coming from, and turns off many useful diagnostics in
  static analysis tools. Exceptions are allowed for modules that purely import
  from an extension to wrap the interface, if it would otherwise cause
  excessive verbosity.
- **Don't create classes with one function**. If your class is an `__init__`
  and a single 'do' function (or a `__call__`), then it can probably be more
  concisely expressed as a function.
- Avoid classes that work by `__call__` unless you have a good reason for it to
  act as a functor. A named function to do the action instead is almost always
  clearer with a proper name.
- **Write docstrings**. We have a mix of styles at the moment, but new
  docstrings should try to follow [Google-style] - it has a good balance
  between clarity and length.
- **Commit messages should be descriptive**; See [How to Write a Git Commit
  Message] for a long explanation of how and why this matters, or just skip to
  [The Seven Rules] for general guidelines (and they are guidelines for us, not
  rules per se). The first line should be a concise summary, ideally not over
  50 characters but never over 72. Follow this first line with an empty line.
  Remember that someone may be looking at your commit in several years time,
  trying to work out the reason for your commit and wondering what on earth you
  were thinking. That someone may be you.

## Pull Requests

- Include a newsfragment before merging. This is a file in the
  `newsfragments` folder with a description of your change. Try to make it
  a one-sentence summary aimed at a DIALS user, not a developer. The file is
  numbered to match the issue or pull request and has 
  [one of the allowed extensions from the list][news-README]. 
  If you don't have an issue number you can also add that file after you've
  created your pull request. Newsfragments are processed on release by a
  tool named [towncrier] and make up our release notes. 
- We aim to squash-merge most pull requests. However, if you are working on a
  longer-term feature branch, or making a lot of changes that might be a
  candidate for a non-squash merge - please try to keep commits a relatively
  clean representation of the implementation of your change, by e.g. using
  `git rebase` - This helps people to understand your changes easier.


## Code linting, automatic formatting and static analysis

- **Please use the pre-commit hooks**. If using the DIALS bootstrap installer
  then this will be done automatically. Otherwise, use `libtbx.precommit
  install` - if in the libtbx ecosystem - or manually install the hooks with
  `pre-commit install`. These use the [pre-commit] tool and ensure that
  various sanity checks are run before commit, including import order,
  formatting, syntax compatibility, basic flake8 checks, lack of conflict
  markers and file size limits. Basically, most of the essential rules will be
  checked automatically by this.
- **We format python code with [black]**. This means that while writing code
  you don't have to worry about laying things out neatly, because black will
  take care of the formatting. We prefer if you commit code formatted with
  black (the pre-commit hook will help do this for you), but if for some reason
  you miss this, the whole codebase is auto-cleaned once a week. Most IDEs
  and editors have support for running formatters like black automatically.
- **Avoid introducing new flake8 warnings** - if you feel that it's appropriate
  to ignore a warning, mark it up explicitly with a [noqa] comment. The most
  important subset of checks are run as part of the pre-commit checks, but
  please try to resolve any other valid warnings shown with a normal run of
  flake8. The configuration in the repository turns off any warnings that
  disagree with our standard practice.
- **We format C++ code with [clang-format]**. We use a configuration for style
  broadly compatible with what our existing prevailing style was. We don't
  require that everyone has clang-format installed - the weekly cleaning job
  will pick it up if you don't - but if you do, remember to run with
  `-style=file` to pick up our configuration.


[pre-commit]: https://github.com/pre-commit/pre-commit
[black]: https://github.com/psf/black
[isort]: https://github.com/PyCQA/isort
[clang-format]: https://clang.llvm.org/docs/ClangFormat.html
[noqa]: http://flake8.pycqa.org/en/3.7.7/user/violations.html#in-line-ignoring-errors
[PEP8]: https://www.python.org/dev/peps/pep-0008
[Google-style]: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
[Zen of Python]: https://www.python.org/dev/peps/pep-0020/#the-zen-of-python
[How to Write a Git Commit Message]: https://chris.beams.io/posts/git-commit
[The Seven Rules]: https://chris.beams.io/posts/git-commit/#seven-rules
[news-README]: https://github.com/dials/dials/blob/master/newsfragments/README.MD
[towncrier]: https://github.com/twisted/towncrier
