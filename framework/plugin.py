class PhilGenerator(object):

    choice_template = '''
    {name}
      .help = "{help}"
    {{
      algorithm = *{choices}
        .type = choice(multi={multi})
        .help = "Select the algorithm to use."
      {parameters}
    }}
    '''

    param_template = '''
      {name}
        .help = "{help}"
      {{
        {parameters}
      }}
    '''

    def __init__(self, mount_point):
        self.mount_point = mount_point

    def __call__(self):

      # Create all the phil parameters
      algorithms = []
      parameters = []
      for name, plugin in self.mount_point.plugins().iteritems():
          algorithms.append(name)
          try:
              parameters.append(self.param_template.format(
                 name = name,
                 help = plugin.__doc__,
                 parameters=  plugin.config))
          except Exception:
              pass

      # Get if multiple choice
      try:
          multi = self.mount_point.multi_choice
      except Exception:
          multi = False

      # Generate the phil string
      phil = self.choice_template.format(
          name = self.mount_point.name,
          help = self.mount_point.__doc__,
          choices = ' '.join(algorithms),
          multi = multi,
          parameters = ' '.join(parameters))

      # Return the phil string
      return phil


class Registry(list):

    def plugins(self):
        plugins = []
        for mp in self:
            for key, value in mp.plugins().iteritems():
                plugins.append(value)
        return plugins

    def configuration(self):
        phil = []
        for mp in self:
            phil.append(mp.configuration())
        return ' '.join(phil)


def localstatic(name, value):
  def decorate(func):
      setattr(func, name, value)
      return func
  return decorate


@localstatic("data", Registry())
def extensions():
    return extensions.data


class Extension(type):

    def __init__(self, name, bases, attrs):
        super(Extension, self).__init__(name, bases, attrs)
        if not hasattr(self, "_plugins"):
            extensions().append(self)
            self._plugins = dict()
        else:
            self._plugins[self.name] = self

    def plugins(self):
        return self._plugins

    def configuration(self):
        return PhilGenerator(self)()

    def factory(self, name, *args, **kwargs):
        return self._plugins[name](*args, **kwargs)
