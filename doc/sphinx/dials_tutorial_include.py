from __future__ import annotations

import os

import sphinx.directives.code as code

"""
Add a directive to literal-include from generated DIALS output logs
and an associated configuration parameter pointing to the location of said logs
"""


def setup(app):
    app.add_config_value("dials_logs", None, "html")
    app.add_directive("dials_tutorial_include", DialsTutorialInclude)
    return {"parallel_read_safe": True}


class DialsTutorialInclude(code.LiteralInclude):
    """Behaves like LiteralInclude but looks for files relative to the known
    DIALS output log directory."""

    def run(self) -> code.List[code.nodes.Node]:
        document = self.state.document
        if not document.settings.file_insertion_enabled:
            return [
                document.reporter.warning("File insertion disabled", line=self.lineno)
            ]
        # convert options['diff'] to absolute path
        if "diff" in self.options:
            _, path = self.env.relfn2path(self.options["diff"])
            self.options["diff"] = path

        try:
            location = self.state_machine.get_source_and_line(self.lineno)
            filename = os.path.join(self.config.dials_logs, self.arguments[0])
            self.env.note_dependency(filename)

            reader = code.LiteralIncludeReader(filename, self.options, self.config)
            text, lines = reader.read(location=location)

            retnode = code.nodes.literal_block(text, text, source=filename)
            self.set_source_info(retnode)
            if self.options.get("diff"):  # if diff is set, set udiff
                retnode["language"] = "udiff"
            elif "language" in self.options:
                retnode["language"] = self.options["language"]
            retnode["linenos"] = (
                "linenos" in self.options
                or "lineno-start" in self.options
                or "lineno-match" in self.options
            )
            retnode["classes"] += self.options.get("class", [])
            extra_args = retnode["highlight_args"] = {}
            if "emphasize-lines" in self.options:
                hl_lines = code.parselinenos(self.options["emphasize-lines"], lines)
                if any(i >= lines for i in hl_lines):
                    code.logger.warning(
                        code.__("line number spec is out of range(1-%d): %r")
                        % (lines, self.options["emphasize-lines"]),
                        location=location,
                    )
                extra_args["hl_lines"] = [x + 1 for x in hl_lines if x < lines]
            extra_args["linenostart"] = reader.lineno_start

            if "caption" in self.options:
                caption = self.options["caption"] or self.arguments[0]
                retnode = code.container_wrapper(self, retnode, caption)

            # retnode will be note_implicit_target that is linked from caption and numref.
            # when options['name'] is provided, it should be primary ID.
            self.add_name(retnode)

            return [retnode]
        except Exception as exc:
            print(exc)
            return [document.reporter.warning(str(exc), line=self.lineno)]
