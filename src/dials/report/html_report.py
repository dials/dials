from __future__ import annotations


class html_report:
    def __init__(self, external_dependencies="remote"):
        self._content = []
        self.external_dependencies = external_dependencies

    def header(self):
        assert self.external_dependencies in ("remote", "local", "embed")

        if self.external_dependencies == "remote":
            plotly_js = (
                '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>'
            )
            bootstrap_js = '<script src="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js"></script>'
            jquery_js = (
                '<script src="https://code.jquery.com/jquery-1.12.0.min.js"></script>'
            )
            bootstrap_css = '<link rel="stylesheet" href="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css">'
            katex_js = '<script type="text/javascript" src=https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.5.1/katex.min.js></script>'
            katex_auto_render_js = '<script type="text/javascript" src=https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.5.1/contrib/auto-render.min.js></script>'
            katex_css = '<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.5.1/katex.min.css">'

        elif self.external_dependencies == "local":
            import libtbx.load_env

            css_dir = libtbx.env.find_in_repositories("dials/static/css")
            js_dir = libtbx.env.find_in_repositories("dials/static/js")
            katex_dir = libtbx.env.find_in_repositories("dials/static/katex")

            plotly_js = f'<script src="{js_dir}/plotly-latest.min.js"></script>'
            bootstrap_js = f'<script src="{js_dir}/bootstrap.min.js"></script>'
            jquery_js = f'<script src="{js_dir}/jquery-1.12.0.min.js"></script>'
            bootstrap_css = (
                f'<link rel="stylesheet" href="{css_dir}/bootstrap.min.css">'
            )
            katex_js = (
                f'<script type="text/javascript" src={katex_dir}/katex.min.js></script>'
            )
            katex_auto_render_js = (
                '<script type="text/javascript" src=%s/contrib/auto-render.min.js></script>'
                % katex_dir
            )
            katex_css = f'<link rel="stylesheet" href="{katex_dir}/katex.min.css">'

        elif self.external_dependencies == "embed":
            import os

            import libtbx.load_env

            css_dir = libtbx.env.find_in_repositories("dials/static/css")
            js_dir = libtbx.env.find_in_repositories("dials/static/js")
            katex_dir = libtbx.env.find_in_repositories("dials/static/katex")

            plotly_js = (
                "<script>%s</script>"
                % open(os.path.join(js_dir, "plotly-latest.min.js")).read()
            )
            bootstrap_js = (
                "<script>%s</script>"
                % open(os.path.join(js_dir, "bootstrap.min.js")).read()
            )
            jquery_js = (
                "<script>%s</script>"
                % open(os.path.join(js_dir, "jquery-1.12.0.min.js")).read()
            )
            bootstrap_css = (
                '<style type="text/css">%s</style>'
                % open(os.path.join(css_dir, "bootstrap.min.css")).read()
            )
            katex_js = (
                "<script>%s</script>"
                % open(os.path.join(katex_dir, "katex.min.js")).read()
            )
            katex_auto_render_js = (
                "<script>%s</script>"
                % open(os.path.join(katex_dir, "contrib/auto-render.min.js")).read()
            )
            katex_css = (
                '<style type="text/css">%s</style>'
                % open(os.path.join(katex_dir, "katex.min.css")).read()
            )

        html_header = f"""
<!DOCTYPE html>
<head>

<!-- Definitions for compatibility with the ccp4i2 browser -->
<script>
    if (!Function.prototype.bind) {{
        Function.prototype.bind = function(oThis) {{
            if (typeof this !== 'function') {{
                // closest thing possible to the ECMAScript 5
                // internal IsCallable function
                throw new TypeError('Function.prototype.bind - what is trying to be bound is not callable');
            }}

            var aArgs   = Array.prototype.slice.call(arguments, 1),
            fToBind = this,
            fNOP    = function() {{}},
            fBound  = function() {{
                return fToBind.apply(this instanceof fNOP
                                     ? this
                                     : oThis,
                                     aArgs.concat(Array.prototype.slice.call(arguments)));
            }};

            if (this.prototype) {{
                // Function.prototype doesn't have a prototype property
                fNOP.prototype = this.prototype;
            }}
            fBound.prototype = new fNOP();

            return fBound;
        }};
    }}
    if (typeof Float64Array === 'undefined') Float64Array = Float32Array;
</script>

<!-- Plotly.js -->
{plotly_js}

<meta name="viewport" content="width=device-width, initial-scale=1" charset="UTF-8">

{jquery_js}
{bootstrap_js}
{katex_js}
{katex_auto_render_js}
{katex_css}
{bootstrap_css}
<style type="text/css">
{self.css()}
</style>

</head>
"""

        return html_header

    def css(self):
        return """
body {
  /*font-family: Helmet, Freesans, Helvetica, Arial, sans-serif;*/
  margin: 8px;
  min-width: 240px;
  margin-left: 5%;
  margin-right: 5%;
}

.plot {
  float: left;
  width: 600px;
  height: 400px;
  margin-bottom: 20px;
}

.katex-display {
  text-align: left;
}

}"""

    def body(self):
        html_body = """

<body>
  %s
  <script>
    renderMathInElement(
        document.body);
  </script>
</body>
""" % ("\n".join(content.html() for content in self._content))

        return html_body

    def html(self):
        return "\n".join((self.header(), self.body()))

    def add_content(self, content):
        self._content.append(content)


class page_header:
    def __init__(self, title):
        self._title = title

    def html(self):
        html = f"""<div class="page-header">
  <h1>{self._title}</h1>
</div>"""
        return html


class panel_group:
    def __init__(self, panels):
        self._panels = panels

    def html(self):
        html = """
<div class="panel-group">
  %s
</div>
""" % ("\n".join(p.html() for p in self._panels))
        return html


class panel:
    def __init__(self, title, panel_id, show=False):
        self._title = title
        self._panel_id = panel_id
        self._content = []
        self._show = show

    def add_content(self, content):
        self._content.append(content)

    def html(self):
        html = """
<div class="panel panel-default">
  <div class="panel-heading" data-toggle="collapse" href="#collapse_%(panel_id)s">
    <h4 class="panel-title">
      <a>%(title)s</a>
    </h4>
  </div>
  <div id="collapse_%(panel_id)s" class="panel-collapse collapse%(in)s">
    <div class="panel-body">
      %(content)s
    </div>
  </div>
</div>
""" % {
            "panel_id": self._panel_id,
            "title": self._title,
            "in": " in" if self._show else "",
            "content": "\n".join(content.html() for content in self._content),
        }
        return html


class table_responsive:
    def __init__(self, table_html, width=None):
        self._table_html = table_html
        self._width = width

    def html(self):
        if self._width is not None:
            style = ' style="width: %ipx"' % self._width
        else:
            style = ""
        html = f"""
<div class="table-responsive"{style}>
  {self._table_html}
</div>
"""
        return html


class plotly_graph:
    def __init__(self, json_data, div_id):
        self._json_data = json_data
        self._div_id = div_id

    def javascript(self):
        import json

        json_str = json.dumps(self._json_data)
        javascript = f"""
  <script>
    var graphs_{self._div_id} = {json_str};
    Plotly.newPlot({self._div_id}, graphs_{self._div_id}.data, graphs_{self._div_id}.layout);
  </script>
  """
        return javascript

    def html(self):
        return "\n".join(
            (
                '<div class="col-xs-6 col-sm-6 col-md-4 plot" id="%s"></div>'
                % self._div_id,
                self.javascript(),
            )
        )


class container_fluid:
    def __init__(self):
        self._content = []

    def add_content(self, content):
        self._content.append(content)

    def html(self):
        html = """
<div class="container-fluid">
%s
</div>
""" % "\n".join(content.html() for content in self._content)
        return html


class div:
    def __init__(self):
        self._content = []

    def add_content(self, content):
        self._content.append(content)

    def html(self):
        html = """
<div>
%s
</div>
""" % "\n".join(content.html() for content in self._content)
        return html


class raw_html:
    def __init__(self, raw_html):
        self._raw_html = raw_html

    def html(self):
        return self._raw_html
