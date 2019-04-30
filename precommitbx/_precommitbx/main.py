import os
import pre_commit.main

# work around https://github.com/Homebrew/homebrew-core/issues/30445
os.environ.pop("__PYVENV_LAUNCHER__", None)


def main():
    return pre_commit.main.main(
        ["run", "--config", ".pre-commit-config.yaml", "--hook-stage", "commit"]
    )


if __name__ == "__main__":
    exit(main())
