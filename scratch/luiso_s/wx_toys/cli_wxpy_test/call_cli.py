import wx
import subprocess

class MyApp(wx.App):
  def OnInit(self):

    print "Hi"
    subprocess.call("dials.python", shell=True)
    print "bye"

    wx.MessageBox("Hi         wx               ", "wxApp")

    return True

if(__name__ == "__main__"):
  app = MyApp(False)
  app.MainLoop()
