#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package."
#
#  learning from an example given by Andrea Gavana
#  This sample shows how to retrieve default platform's
#  bitmaps using wx.ArtProvider

import wx

class BitmapFrame(wx.Frame):

    def __init__(self):

        wx.Frame.__init__(self, None, -1, title='ArtProvider example')

        panel = wx.Panel(self)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        bitmap_sizer = wx.BoxSizer(wx.HORIZONTAL)

        bitmap_sizer.Add((0, 0), 1, wx.EXPAND)

        # Show a few bitmaps retrieved via wx.ArtProvider
        for kind in [wx.ART_GO_UP, wx.ART_GO_DOWN \
                     , wx.ART_GO_BACK
                     , wx.ART_GO_FORWARD
                     , wx.ART_HELP_PAGE
                     , wx.ART_HELP_BOOK
                     , wx.ART_HELP_FOLDER
                     , wx.ART_FILE_OPEN
                     , wx.ART_FILE_SAVE
                     , wx.ART_FILE_SAVE_AS
                     , wx.ART_PRINT
                     , wx.ART_HELP
                     , wx.ART_TIP
                     , wx.ART_REPORT_VIEW
                     , wx.ART_LIST_VIEW
                     , wx.ART_NEW_DIR
                     , wx.ART_HARDDISK
                     , wx.ART_FLOPPY
                     , wx.ART_CDROM
                     , wx.ART_REMOVABLE
                     , wx.ART_FOLDER
                     , wx.ART_FOLDER_OPEN
                     , wx.ART_GO_DIR_UP
                     , wx.ART_EXECUTABLE_FILE
                     , wx.ART_NORMAL_FILE
                     , wx.ART_TICK_MARK
                     , wx.ART_CROSS_MARK
                     , wx.ART_ERROR
                     , wx.ART_QUESTION
                     , wx.ART_WARNING
                     , wx.ART_INFORMATION
                     , wx.ART_MISSING_IMAGE
                     , wx.ART_COPY
                     , wx.ART_CUT
                     , wx.ART_PASTE
                     , wx.ART_DELETE
                     , wx.ART_NEW
                     , wx.ART_UNDO
                     , wx.ART_REDO
                     , wx.ART_QUIT
                     , wx.ART_FIND
                     , wx.ART_FIND_AND_REPLACE
                     ]:
            bmp = wx.ArtProvider.GetBitmap(kind, wx.ART_OTHER, (32, 32))
            static_bitmap = wx.StaticBitmap(panel, -1, bmp)
            bitmap_sizer.Add(static_bitmap, 0, wx.ALL, 5)

        # Layout everything in a nice sizer
        bitmap_sizer.Add((0, 0), 1, wx.EXPAND)
        main_sizer.Add((0, 0), 1, wx.EXPAND)
        main_sizer.Add(bitmap_sizer, 0, wx.EXPAND)
        main_sizer.Add((0, 0), 1, wx.EXPAND)

        panel.SetSizer(main_sizer)
        main_sizer.SetSizeHints(panel)
        main_sizer.Layout()


if __name__ == '__main__':
    app = wx.App(0)
    frame = BitmapFrame()
    frame.Show()
    app.MainLoop()