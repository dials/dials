towncrier has a few standard types of news fragments,
signified by the file extension. These are:

    .feature: Signifying a new feature.
    .bugfix: Signifying a bug fix.
    .doc: Signifying a documentation improvement.
    .removal: Signifying a deprecation or removal of public API.
    .misc: A ticket has been closed, but it is not of interest to users.

The start of the filename is the ticket number, and the content is what will
end up in the news file. For example, if ticket #850 is about adding a new
widget, the filename would be `newsfragments/850.feature` and the
content would be
```
myproject.widget has been added.
```
