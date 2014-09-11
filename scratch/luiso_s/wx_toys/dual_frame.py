import wx

class ChildFrame(wx.Frame):
    def __init__(self, parent):
        wx.Frame.__init__(self, None, size=(150,100), title='ChildFrame')
        self.parent = parent
        pan = wx.Panel(self)
        self.txt = wx.TextCtrl(pan, -1, pos=(0,0), size=(500,20), style=wx.DEFAULT)
        self.but = wx.Button(pan,-1, pos=(10,30), label='Tell parent')
        self.Bind(wx.EVT_BUTTON, self.onbutton, self.but)

    def onbutton(self, evt):
        text = self.txt.GetValue()
        self.parent.txt.write('Child says: %s' %text)

class MainFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, size=(150,100), title='MainFrame')
        pan =wx.Panel(self)
        self.txt = wx.TextCtrl(pan, -1, pos=(0,0), size=(500,20), style=wx.DEFAULT)
        self.but = wx.Button(pan,-1, pos=(10,30), label='Tell child')
        self.Bind(wx.EVT_BUTTON, self.onbutton, self.but)
        self.child = ChildFrame(self)
        self.child.Show()

    def onbutton(self, evt):
        text = self.txt.GetValue()
        self.child.txt.write('Parent says: %s' %text)





if __name__ == "__main__":

    App=wx.PySimpleApp()
    MainFrame().Show()
    App.MainLoop()

