import subprocess
import wx


class element(object):
  pass


class MyWidg(wx.Frame):
  def __init__(self, parent, id, title):
    super(MyWidg, self).__init__(parent, id, title)

    wx_btn_lst = []
    self.step_btn_lst = []

    import_btn = element()
    import_btn.lbl = "\n import \n"
    import_btn.cmd = "dials.import"
    self.step_btn_lst.append(import_btn)

    find_btn = element()
    find_btn.lbl = "\n find_spots \n"
    find_btn.cmd = "dials.find_spots"
    self.step_btn_lst.append(find_btn)


    index_btn = element()
    index_btn.lbl = "\n index \n"
    index_btn.cmd = "dials.index"
    self.step_btn_lst.append(index_btn)


    for step_btn in self.step_btn_lst:
      exmpl_btt = wx.Button(self, -1, step_btn.lbl, (0,0))
      wx_btn_lst.append(exmpl_btt)
      step_btn.bt_id = exmpl_btt.GetId()

    main_panel = wx.Panel(self)

    go_btn = wx.Button(self, -1, " Go ", (0,0))
    self.cli_txt = wx.TextCtrl(self, -1, size = (800,20), style = wx.DEFAULT)


    vbox = wx.BoxSizer(wx.VERTICAL)

    for wx_btn in wx_btn_lst:
      vbox.Add(wx_btn, 0, border = 3)


    hbox = wx.BoxSizer(wx.HORIZONTAL)
    hbox.Add(vbox, 0, border = 3)
    hbox.Add(main_panel, 0, border = 3)

    bt_box = wx.BoxSizer(wx.HORIZONTAL)
    bt_box.Add(self.cli_txt,wx.EXPAND | wx.ALL, border = 3)
    bt_box.Add(go_btn, 0, border = 3)

    bg_box = wx.BoxSizer(wx.VERTICAL)
    bg_box.Add(hbox, 0, border = 3)
    bg_box.Add(bt_box, 0, border = 3)

    self.SetSizer(bg_box)
    bg_box.Fit(self)

    for wx_btn in wx_btn_lst:
      wx_btn.Bind(wx.EVT_BUTTON, self.on_btn)

    go_btn.Bind(wx.EVT_BUTTON, self.on_go_button)


  def on_btn(self, event):

    for step_btn in self.step_btn_lst:
      if(step_btn.bt_id == event.GetId()):
        self.cli_txt.SetValue(step_btn.cmd)


  def on_go_button(self, event):
    txt_to_run = self.cli_txt.GetValue()
    subprocess.call(txt_to_run, shell=True)


class MyApp(wx.App):
  def OnInit(self):


    frame = MyWidg(None, -1, 'menu1.py')
    frame.Show(True)
    return True

app = MyApp(0)
app.MainLoop()
